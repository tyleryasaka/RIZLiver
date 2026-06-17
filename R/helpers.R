.estimate_nb_size = function(counts, min_mean = 0.5, q = 0.5, fallback = 10) {
  X = as.matrix(counts)
  if (nrow(X) < 50) { warning('too few genes to estimate dispersion; using fallback r0'); return(fallback) }
  n = ncol(X)
  N = colSums(X); N[N == 0] = 1
  cvN2 = var(N) / mean(N)^2
  mu = rowMeans(X)
  v  = rowSums((X - mu)^2) / (n - 1)
  keep = which(mu >= min_mean & v > mu)
  if (length(keep) < 20) { warning('too few usable genes to estimate dispersion; using fallback r0'); return(fallback) }
  phi_naive = (v[keep] - mu[keep]) / mu[keep]^2
  phi = (quantile(phi_naive, q, names = FALSE) - cvN2) / (1 + cvN2)
  if (!is.finite(phi) || phi <= 0) { warning('non-positive dispersion estimate; using fallback r0'); return(fallback) }
  1 / phi
}

.build_reference_for_fitting = function(refs, interp_points = 100) {
  ref_per_source = list()
  positions = NULL

  for (src_idx in seq_along(refs)) {
    m = refs[[src_idx]]
    src_name = if (!is.null(names(refs)[src_idx])) names(refs)[src_idx] else paste0('s', src_idx)

    if (is.data.frame(m)) m = as.matrix(m)
    if (is.null(dim(m)) || nrow(m) < 10 || ncol(m) < 2) {
      warning(sprintf("source '%s': not a genes x regions matrix with >=2 regions; skipping", src_name))
      next
    }
    m = m[!apply(m, 1, anyNA), , drop = FALSE]
    if (nrow(m) < 10) {
      warning(sprintf("source '%s': fewer than 10 genes after NA removal; skipping", src_name))
      next
    }

    n_reg = ncol(m)
    region_pos = (seq_len(n_reg) - 1) / (n_reg - 1)
    xout = seq(0, 1, length.out = interp_points)
    src_interp = t(apply(m, 1, function(x) approx(region_pos, x, xout = xout)$y))
    rownames(src_interp) = rownames(m)
    colnames(src_interp) = paste0(src_name, '_p', seq_len(interp_points))

    ref_per_source[[src_name]] = src_interp
    if (is.null(positions)) positions = xout
  }

  if (length(ref_per_source) < 1) stop('no usable reference sources')

  all_genes = Reduce(union, lapply(ref_per_source, rownames))
  if (length(all_genes) < 10)
    stop(sprintf('only %d union genes across sources', length(all_genes)))

  ref_per_source = lapply(ref_per_source, function(m) {
    out = matrix(NA_real_, nrow = length(all_genes), ncol = ncol(m),
                 dimnames = list(all_genes, colnames(m)))
    out[rownames(m), ] = m
    out
  })

  list(refs = ref_per_source, positions = positions)
}

.fit_layer_profiles = function(refs_per_source, positions, landmark_genes, r0,
                               n_layers = 9, layer_prior = NULL) {
  if (length(refs_per_source) < 1) stop('no reference sources')
  stopifnot(length(positions) == ncol(refs_per_source[[1]]))

  lm = intersect(landmark_genes, rownames(refs_per_source[[1]]))
  if (length(lm) < 2) stop('fewer than 2 landmark genes present in the reference')

  br = seq(min(positions), max(positions), length.out = n_layers + 1)
  layer = cut(positions, breaks = br, labels = FALSE, include.lowest = TRUE)
  occ = tabulate(layer, nbins = n_layers)

  n_sources = length(refs_per_source)
  per_src_lm = array(NA_real_, dim = c(length(lm), n_layers, n_sources),
                     dimnames = list(lm, NULL, names(refs_per_source)))

  for (s in seq_len(n_sources)) {
    M = refs_per_source[[s]]
    Mtot = Matrix::colSums(M, na.rm = TRUE); Mtot[Mtot == 0] = 1
    Mfrac = sweep(as.matrix(M[lm, , drop = FALSE]), 2, Mtot, '/')
    for (k in seq_len(n_layers)) {
      ck = which(layer == k)
      if (length(ck) < 1) next
      per_src_lm[, k, s] = rowMeans(Mfrac[, ck, drop = FALSE], na.rm = TRUE)
    }
  }

  alpha = beta = matrix(NA_real_, length(lm), n_layers, dimnames = list(lm, NULL))
  for (k in seq_len(n_layers)) {
    for (g_idx in seq_len(length(lm))) {
      v = per_src_lm[g_idx, k, ]
      v = v[is.finite(v) & v > 0]
      n_s = length(v)
      if (n_s < 1) {
        alpha[g_idx, k] = r0
        beta[g_idx, k] = r0 / 1e-10
        next
      }
      m_hat = max(mean(v), 1e-12)
      if (n_s >= 2) {
        sigma2 = var(v)
        if (is.finite(sigma2) && sigma2 > 1e-20) {
          alpha[g_idx, k] = m_hat^2 / sigma2
          beta[g_idx, k]  = m_hat / sigma2
          next
        }
      }
      alpha[g_idx, k] = r0
      beta[g_idx, k]  = r0 / m_hat
    }
  }

  ok = which(occ > 0)
  if (length(ok) < 2) stop('fewer than 2 usable layers; check reference / n_layers')
  prior = if (is.null(layer_prior)) occ / sum(occ) else layer_prior / sum(layer_prior)
  list(lm = lm, alpha = alpha, beta = beta, ok = ok,
       logprior = log(prior + 1e-12), n_layers = n_layers, r0 = r0)
}

.update_layer_profiles_bayes = function(counts, zonation, prior_fit, n_layers) {
  bin = pmin(pmax(ceiling(zonation * n_layers), 1), n_layers)
  N = Matrix::colSums(counts); N[N == 0] = 1

  lm = rownames(prior_fit$alpha)
  lm_in_counts = intersect(lm, rownames(counts))
  if (length(lm_in_counts) < 2) stop('fewer than 2 landmark genes available for second pass')

  X = as.matrix(counts[lm_in_counts, , drop = FALSE])
  idx = match(lm_in_counts, lm)

  alpha = prior_fit$alpha
  beta  = prior_fit$beta

  for (k in seq_len(n_layers)) {
    cells_k = which(bin == k)
    if (length(cells_k) < 1) next
    sum_x = rowSums(X[, cells_k, drop = FALSE])
    sum_N = sum(N[cells_k])
    alpha[idx, k] = alpha[idx, k] + sum_x
    beta[idx, k]  = beta[idx, k]  + sum_N
  }

  occ = tabulate(bin, nbins = n_layers)
  ok = which(occ > 0)
  if (length(ok) < 2) stop('fewer than 2 usable layers in second pass')
  prior = occ / sum(occ)
  list(lm = lm, alpha = alpha, beta = beta, ok = ok,
       logprior = log(prior + 1e-12), n_layers = n_layers, r0 = prior_fit$r0)
}

.expected_layer = function(counts, fit) {
  g = intersect(fit$lm, rownames(counts))
  if (length(g) < 2) {
    warning('fewer than 2 landmark genes present in counts')
    return(setNames(rep(NA_real_, ncol(counts)), colnames(counts)))
  }
  alpha = fit$alpha[g, , drop = FALSE]
  beta  = fit$beta[g, , drop = FALSE]
  r0 = fit$r0
  ok = fit$ok
  logprior = fit$logprior
  N = Matrix::colSums(counts); N[N == 0] = 1
  X = as.matrix(counts[g, , drop = FALSE])

  el = vapply(seq_len(ncol(counts)), function(j) {
    ll = vapply(ok, function(k) {
      a = alpha[, k]; b = beta[, k]
      size_eff = a * r0 / (a + r0)
      prob_eff = b * r0 / (b * r0 + N[j] * (a + r0))
      sum(dnbinom(X[, j], size = size_eff, prob = prob_eff, log = TRUE))
    }, numeric(1))
    ll = ll + logprior[ok]; ll = ll - max(ll)
    post = exp(ll); sum(ok * post / sum(post))
  }, numeric(1))
  setNames(el, colnames(counts))
}

.checkInputs = function(counts, refs_per_source, positions, landmark_genes) {
  if (is.null(rownames(counts)) || is.null(colnames(counts)))
    stop('counts must have gene names as rownames and cell/spot names as colnames.')
  if (min(counts, na.rm = TRUE) < 0)
    stop('counts contains negative values; raw non-negative counts are expected.')
  if (length(refs_per_source) < 1)
    stop('refs_per_source must contain at least one source.')
  ref_genes = Reduce(union, lapply(refs_per_source, rownames))
  if (is.null(ref_genes))
    stop('reference sources must have gene names as rownames.')
  if (length(positions) != ncol(refs_per_source[[1]]))
    stop('positions must have one value per column of each reference source.')
  lm = Reduce(intersect, list(landmark_genes, ref_genes, rownames(counts)))
  if (length(lm) < 2)
    stop('fewer than 2 landmark genes shared across landmark_genes, reference, and counts; check gene naming/species.')
  if (length(lm) < 4)
    warning(sprintf('only %d landmark genes usable; estimate may be weak.', length(lm)))
  invisible(NULL)
}

fitZonation = function(counts, refs_per_source, positions, landmark_genes,
                       n_layers = 9, layer_prior = NULL,
                       r0 = NULL, disp_min_mean = 0.5, disp_quantile = 0.5,
                       second_pass = TRUE, verbose = FALSE) {
  .checkInputs(counts, refs_per_source, positions, landmark_genes)
  cells = colnames(counts)

  if (is.null(r0)) {
    r0 = .estimate_nb_size(counts, min_mean = disp_min_mean, q = disp_quantile)
    if (verbose) cat(sprintf('estimated NB size r0 = %.4g\n', r0))
  }

  fit_1 = .fit_layer_profiles(refs_per_source, positions, landmark_genes, r0, n_layers, layer_prior)
  s_base_1 = .expected_layer(counts, fit_1)
  cal_1 = approxfun(x = s_base_1, y = dplyr::percent_rank(s_base_1), rule = 2)
  label_1 = setNames(cal_1(s_base_1), cells)
  if (verbose) cat(sprintf('first pass: %d landmark genes, %d layers\n',
                           length(fit_1$lm), length(fit_1$ok)))

  if (!second_pass) return(list(fit = fit_1, cal = cal_1, label = label_1))

  fit_2 = .update_layer_profiles_bayes(counts, label_1, fit_1, n_layers)
  s_base_2 = .expected_layer(counts, fit_2)
  cal_2 = approxfun(x = s_base_2, y = dplyr::percent_rank(s_base_2), rule = 2)
  label_2 = setNames(cal_2(s_base_2), cells)
  if (verbose)
    cat(sprintf('second pass: %d landmark genes updated, %d layers\n',
                length(fit_2$lm), length(fit_2$ok)))

  list(fit = fit_2, cal = cal_2, label = label_2)
}

predictZonation = function(model, counts) {
  cells = colnames(counts)
  s = .expected_layer(counts, model$fit)
  setNames(model$cal(s), cells)
}
