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

.expected_layer = function(counts, fit) {
  g = intersect(fit$lm, rownames(counts))
  if (length(g) < 2) {
    warning('fewer than 2 landmark genes present in counts')
    return(setNames(rep(NA_real_, ncol(counts)), colnames(counts)))
  }
  alpha = fit$alpha[g, , drop = FALSE]; beta = fit$beta[g, , drop = FALSE]
  ok = fit$ok; logprior = fit$logprior
  N = Matrix::colSums(counts); N[N == 0] = 1
  X = as.matrix(counts[g, , drop = FALSE])
  el = vapply(seq_len(ncol(counts)), function(j) {
    ll = vapply(ok, function(k)
      sum(dnbinom(X[, j], size = alpha[, k], prob = beta[, k] / (beta[, k] + N[j]), log = TRUE)),
      numeric(1))
    ll = ll + logprior[ok]; ll = ll - max(ll)
    post = exp(ll); sum(ok * post / sum(post))
  }, numeric(1))
  setNames(el, colnames(counts))
}

.checkInputs = function(counts, reference, reference_pos, landmark_genes) {
  if (is.null(rownames(counts)) || is.null(colnames(counts)))
    stop('counts must have gene names as rownames and cell/spot names as colnames.')
  if (min(counts, na.rm = TRUE) < 0)
    stop('counts contains negative values; raw non-negative counts are expected.')
  if (is.null(rownames(reference)))
    stop('reference must have gene names as rownames.')
  if (length(reference_pos) != ncol(reference))
    stop('reference_pos must have one value per column of reference.')
  lm = Reduce(intersect, list(landmark_genes, rownames(reference), rownames(counts)))
  if (length(lm) < 2)
    stop('fewer than 2 landmark genes shared across landmark_genes, reference, and counts; check gene naming/species.')
  if (length(lm) < 4)
    warning(sprintf('only %d landmark genes usable; estimate may be weak.', length(lm)))
  invisible(NULL)
}

.update_layer_profiles_bayes = function(counts, zonation, discovered_genes, prior_fit,
                                        r0, w_prior = 10, n_layers) {
  bin = pmin(pmax(ceiling(zonation * n_layers), 1), n_layers)
  N = Matrix::colSums(counts); N[N == 0] = 1

  lm_with_prior = intersect(rownames(prior_fit$alpha), rownames(counts))
  all_g = unique(c(lm_with_prior, discovered_genes))
  if (length(all_g) < 2) stop('fewer than 2 genes available for second pass')

  X = as.matrix(counts[all_g, , drop = FALSE])
  has_prior = all_g %in% lm_with_prior

  alpha = beta = matrix(NA_real_, length(all_g), n_layers,
                        dimnames = list(all_g, NULL))
  occ = numeric(n_layers)

  for (k in seq_len(n_layers)) {
    cells_k = which(bin == k); occ[k] = length(cells_k)
    if (length(cells_k) < 1) next
    sum_x = rowSums(X[, cells_k, drop = FALSE])
    sum_N = sum(N[cells_k])

    # Default posterior mean: data-only (no informative prior)
    m_post = (sum_x + 1e-12) / (sum_N + 1e-12)

    # Conjugate Gamma-Poisson update for genes that have a reference prior
    if (any(has_prior)) {
      g_p   = all_g[has_prior]
      m_ref = prior_fit$alpha[g_p, k] / prior_fit$beta[g_p, k]
      valid = !is.na(m_ref) & m_ref > 0
      if (any(valid)) {
        idx     = which(has_prior)[valid]
        m_ref_v = m_ref[valid]
        a_post  = w_prior * r0 + sum_x[idx]
        b_post  = w_prior * r0 / m_ref_v + sum_N
        m_post[idx] = a_post / b_post
      }
    }

    m_post[m_post < 1e-12] = 1e-12
    alpha[, k] = r0
    beta[, k]  = r0 / m_post
  }

  ok = which(!is.na(alpha[1, ]))
  if (length(ok) < 2) stop('fewer than 2 usable layers in second pass')

  prior = occ / sum(occ)
  list(lm = all_g, alpha = alpha, beta = beta, ok = ok,
       logprior = log(prior + 1e-12), n_layers = n_layers, r0 = r0)
}

.build_reference_for_fitting = function(refs, interp_points = 100) {
  ref_cols    = list()
  ref_pos_all = numeric(0)

  for (src_idx in seq_along(refs)) {
    m        = refs[[src_idx]]
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

    n_reg      = ncol(m)
    region_pos = (seq_len(n_reg) - 1) / (n_reg - 1)
    xout       = seq(0, 1, length.out = interp_points)
    src_interp = t(apply(m, 1, function(x) approx(region_pos, x, xout = xout)$y))
    rownames(src_interp) = rownames(m)
    colnames(src_interp) = paste0(src_name, '_p', seq_len(interp_points))

    ref_cols[[length(ref_cols) + 1]] = src_interp
    ref_pos_all = c(ref_pos_all, xout)
  }

  if (length(ref_cols) < 1) stop('no usable reference sources')

  all_genes = Reduce(union, lapply(ref_cols, rownames))
  if (length(all_genes) < 10)
    stop(sprintf('only %d union genes across sources', length(all_genes)))

  ref_cols = lapply(ref_cols, function(m) {
    out = matrix(NA_real_, nrow = length(all_genes), ncol = ncol(m),
                 dimnames = list(all_genes, colnames(m)))
    out[rownames(m), ] = m
    out
  })

  list(matrix = do.call(cbind, ref_cols), positions = ref_pos_all)
}

.fit_layer_profiles = function(reference, reference_pos, landmark_genes, r0,
                               n_layers = 9, layer_prior = NULL) {
  stopifnot(length(reference_pos) == ncol(reference))
  lm = intersect(landmark_genes, rownames(reference))
  if (length(lm) < 2) stop('fewer than 2 landmark genes present in the reference')
  br = seq(min(reference_pos), max(reference_pos), length.out = n_layers + 1)
  layer = cut(reference_pos, breaks = br, labels = FALSE, include.lowest = TRUE)
  Rtot = Matrix::colSums(reference, na.rm = TRUE); Rtot[Rtot == 0] = 1
  Rfrac = sweep(as.matrix(reference[lm, , drop = FALSE]), 2, Rtot, "/")
  alpha = beta = matrix(NA_real_, length(lm), n_layers, dimnames = list(lm, NULL))
  occ = numeric(n_layers)
  for (k in seq_len(n_layers)) {
    ck = which(layer == k); occ[k] = length(ck)
    if (length(ck) < 1) next
    m = rowMeans(Rfrac[, ck, drop = FALSE], na.rm = TRUE)
    m[!is.finite(m) | m < 1e-12] = 1e-12
    alpha[, k] = r0
    beta[, k]  = r0 / m
  }
  ok = which(!is.na(alpha[1, ]))
  if (length(ok) < 2) stop('fewer than 2 usable layers; check reference_pos / n_layers')
  prior = if (is.null(layer_prior)) occ / sum(occ) else layer_prior / sum(layer_prior)
  list(lm = lm, alpha = alpha, beta = beta, ok = ok,
       logprior = log(prior + 1e-12), n_layers = n_layers, r0 = r0)
}

fitZonation = function(counts, reference, reference_pos, landmark_genes,
                       n_layers = 9, layer_prior = NULL,
                       r0 = NULL, disp_min_mean = 0.5, disp_quantile = 0.5,
                       second_pass = TRUE, w_prior = 10, verbose = FALSE) {
  .checkInputs(counts, reference, reference_pos, landmark_genes)
  cells = colnames(counts)

  if (is.null(r0)) {
    r0 = .estimate_nb_size(counts, min_mean = disp_min_mean, q = disp_quantile)
    if (verbose) cat(sprintf('estimated NB size r0 = %.4g\n', r0))
  }

  fit_1 = .fit_layer_profiles(reference, reference_pos, landmark_genes, r0, n_layers, layer_prior)
  s_base_1 = .expected_layer(counts, fit_1)
  cal_1 = approxfun(x = s_base_1, y = dplyr::percent_rank(s_base_1), rule = 2)
  label_1 = setNames(cal_1(s_base_1), cells)
  if (verbose) cat(sprintf('first pass: %d landmark genes, %d layers\n',
                           length(fit_1$lm), length(fit_1$ok)))

  if (!second_pass) return(list(fit = fit_1, cal = cal_1, label = label_1))

  fit_2 = .update_layer_profiles_bayes(counts, label_1, character(0), fit_1,
                                       r0, w_prior, n_layers)
  s_base_2 = .expected_layer(counts, fit_2)
  cal_2 = approxfun(x = s_base_2, y = dplyr::percent_rank(s_base_2), rule = 2)
  label_2 = setNames(cal_2(s_base_2), cells)
  if (verbose)
    cat(sprintf('second pass: %d landmark genes updated (w_prior=%.2f), %d layers\n',
                length(fit_2$lm), w_prior, length(fit_2$ok)))

  list(fit = fit_2, cal = cal_2, label = label_2)
}

predictZonation = function(model, counts) {
  cells = colnames(counts)
  s = .expected_layer(counts, model$fit)
  setNames(model$cal(s), cells)
}
