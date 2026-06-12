ZONES = c('Zone_1', 'Zone_2', 'Zone_3')
ZONE_COLORS = c('#504880FF', '#1BB6AFFF', '#FFAD0AFF', '#D72000FF', '#F080D8FF')
minMaxNorm = function(v) (v - min(v, na.rm = T)) / (max(v, na.rm = T) - min(v, na.rm = T))
.RIZLiver_cache = new.env(parent = emptyenv())
.obj_id_counter = new.env(parent = emptyenv())
.obj_id_counter$n = 0L

next_obj_id = function() {
  .obj_id_counter$n <- .obj_id_counter$n + 1L
  sprintf("%.6f_%d", as.numeric(Sys.time()), .obj_id_counter$n)
}

access_cache = function(zone_obj, prop, fn) {
  id = zone_obj$obj_id
  if (is.null(.RIZLiver_cache[[id]])) {
    .RIZLiver_cache[[id]] = new.env(parent = emptyenv())
  }
  if (is.null(.RIZLiver_cache[[id]][[prop]])) {
    .RIZLiver_cache[[id]][[prop]] = fn()
  }
  .RIZLiver_cache[[id]][[prop]]
}

cleanMatrix = function(mtx) {
  mtx = as.matrix(mtx)
  mtx[is.na(mtx)] = 0
  mtx
}

predict_position = function(zone_obj) {
  access_cache(zone_obj, 'predict_position', function() {
    predictZonation(zone_obj, zone_obj$mtx)
  })
}

getZonationGradient_help = function(zone_obj) {
  predict_position(zone_obj) * 2 + 1
}

#' Obtain an interpolation data frame
#'
#' @param coords Coordinate matrix with samples as rows, and columns `x` and `y`. Rownames of coords should match colnames of zone_obj$mtx.
#' @param zone_obj Calibrated Zonation Object
#' @param resolution Optional numeric value for the resolution, where higher value results in a more granular interpolation (default 1)
#' @return An interpolation data frame
#' @noRd
apply_interpolation = function(coords, zone_obj, resolution = 1) {
  coords$zonation = getZonationGradient_help(zone_obj)[rownames(coords)]
  nx = abs((range(coords$x)[[1]] - range(coords$x)[[2]]) / 300 * resolution)
  ny = abs((range(coords$y)[[1]] - range(coords$y)[[2]]) / 300 * resolution)
  with(coords, akima::interp(x, y, zonation, duplicate = "mean", linear = TRUE, extrap = FALSE, nx = nx, ny = ny))
}

#' Train a zonation model on a baseline sample
#'
#' @param mtx Gene expression matrix (raw counts) with genes as rows.
#' @param coords (Optional) Coordinate matrix with samples as rows and columns `x` and `y`. Used to set the spatial scale factor.
#' @param species (Optional) Species to use, defaults to `'human'`. Supports `'human'` and `'mouse'`.
#' @param eval (Optional) If TRUE, only use the liver cell atlas reference for that species, excluding other sources. Default FALSE.
#' @param verbose (Optional) Print fit diagnostics.
#' @return A trained zonation object: a list with `fit`, `cal`, `label`, `mtx`, `obj_id`, `scale_factor`.
#' @export
trainModel = function(mtx, coords = NULL, species = 'human', eval = FALSE, verbose = FALSE) {
  if (species == 'human') {
    initial_weights = readRDS(system.file('extdata', 'initial_weights_human.RDS', package = 'RIZLiver'))
  } else if (species == 'mouse') {
    initial_weights = readRDS(system.file('extdata', 'initial_weights_mouse.RDS', package = 'RIZLiver'))
  } else {
    stop("Only 'human' and 'mouse' species are supported at the moment. (Specify with species = 'mouse'")
  }
  initial_weights = initial_weights[abs(initial_weights) > 0.2]
  initial_weights_zone_1 = initial_weights[initial_weights < 0]
  initial_weights_zone_3 = initial_weights[initial_weights > 0]
  max_markers = 100
  zone_1.landmark = names(head(initial_weights_zone_1[order(initial_weights_zone_1)], max_markers))
  zone_3.landmark = names(tail(initial_weights_zone_3[order(initial_weights_zone_3)], max_markers))

  ref_dir = system.file('extdata', 'references', species, package = 'RIZLiver')
  if (!nzchar(ref_dir) || !dir.exists(ref_dir))
    stop(sprintf("No reference directory found for species '%s'", species))

  if (eval) {
    eval_file = file.path(ref_dir, sprintf('liver_cell_%s.RDS', species))
    if (!file.exists(eval_file))
      stop(sprintf("eval=TRUE but %s not found", eval_file))
    ref_files = eval_file
  } else {
    ref_files = list.files(ref_dir, pattern = '\\.RDS$', full.names = TRUE, ignore.case = TRUE)
    if (length(ref_files) < 1)
      stop(sprintf("No reference .RDS files found in %s", ref_dir))
  }
  refs = setNames(lapply(ref_files, readRDS),
                  tools::file_path_sans_ext(basename(ref_files)))
  if (verbose)
    cat(sprintf('loaded %d reference source(s) for %s: %s\n',
                length(refs), species, paste(names(refs), collapse = ', ')))

  mtx = cleanMatrix(mtx)
  combined = .build_reference_for_fitting(refs)

  zone_obj = fitZonation(mtx, combined$matrix, combined$positions,
                         c(zone_1.landmark, zone_3.landmark), verbose = verbose)
  applyModel(mtx, zone_obj, coords = coords)
}

#' Apply a trained zonation model to a sample
#'
#' Attaches a count matrix (and optional coordinates) to a fitted zonation object so
#' downstream functions can use it. Called internally by `trainModel`; call directly
#' to apply an already-fitted model to a new sample.
#'
#' @param mtx Gene expression matrix with genes as rows and samples as columns
#' @param zone_obj A fitted zonation object (output of `fitZonation` or `trainModel`)
#' @param coords (Optional) Coordinate matrix with samples as rows and columns `x` and `y`. Used to set the spatial scale factor for downstream interpolation.
#' @return The zonation object with `mtx`, `obj_id`, and `scale_factor` set.
#' @export
applyModel = function(mtx, zone_obj, coords = NULL) {
  mtx = cleanMatrix(mtx)
  zone_obj$mtx = mtx
  zone_obj$obj_id = next_obj_id()
  zone_obj$scale_factor = if (is.null(coords)) 1 else {
    x_range = abs(diff(range(coords$x)))
    if (x_range < 100) 50 else 1
  }
  zone_obj
}

#' Per-layer mean expression fraction for each landmark gene, learned from the reference
#'
#' @param zone_obj Calibrated zonation object
#' @return A data frame with one row per landmark gene and one column per layer.
#'   Values are the mean expression fractions used in the model.
#' @export
getGeneZonation = function(zone_obj) {
  m = zone_obj$fit$r0 / zone_obj$fit$beta
  rownames(m) = zone_obj$fit$lm
  colnames(m) = paste0('layer_', seq_len(ncol(m)))
  as.data.frame(m)
}

#' Apply the model to new samples, returning the zonation per cell/spot as a gradient between 1 and 3
#'
#' @param zone_obj Calibrated Zonation Object
#' @return A vector of numeric zonation assignments (continuous)
#' @export
getZonationGradient = function(zone_obj) {
  getZonationGradient_help(zone_obj)
}

#' Apply the model to new samples, returning the zonation per cell/spot as discrete bins (1, 2, or 3)
#'
#' @param zone_obj Calibrated Zonation Object
#' @return A vector of zonation assignments (discrete)
#' @export
getZone = function(zone_obj) {
  mtx = zone_obj$mtx
  zonescore = predict_position(zone_obj)
  zone = ifelse(zonescore < (1 / 3), 'Zone_1', ifelse(zonescore < (2 / 3), 'Zone_2', 'Zone_3'))
  zone = factor(zone, levels = ZONES)
  names(zone) = colnames(mtx)
  zone
}

#' Plot a gene's expression along the zonation axis with a polynomial regression line
#'
#' @param zone_obj Calibrated zonation object
#' @param gene Name of the gene to plot
#' @return A ggplot object
#' @export
plotRegression = function(zone_obj, gene) {
  mtx = zone_obj$mtx
  df = data.frame(zonation = getZonationGradient(zone_obj))
  df[[gene]] = mtx[gene, ]
  ggplot(df, aes(x = zonation, y = .data[[gene]])) +
    geom_point(color = '#504880FF') +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = TRUE, color = "#F080D8FF", fill = "#C8C0F8FF") +
    xlab('Zonation')
}

#' Use this function to infer zonation of all cell types based on their spatial interpolation within the zonation gradient.
#' Zonation gradient is inferred from a subset of samples (if specified; e.g. hepatocytes) and then applied to all samples via interpolation.
#' This function will return an output very similar to getZone, but may be slightly "smoothened".
#' An additional benefit is that this can be used to infer zonation from hepatocytes, then applied to non-parenchymal cells based on proximity to zonated hepatocytes, rather than their own gene expresison.
#'
#' @param coords Coordinate matrix with samples as rows, and columns `x` and `y`. Rownames of coords should match colnames of zone_obj$mtx.
#' @param zone_obj Calibrated Zonation Object
#' @param resolution Optional numeric value for the resolution, where higher value results in a more granular interpolation (default 1)
#' @param use_for_inference (optional) A vector of sample names which should be used for zonation inference (recommended to use only hepatocytes, if annotation is available). If not provided, all samples will be used.
#' @return A vector of zonation assignments (discrete) for all samples
#' @export
getZoneSpatial = function(coords, zone_obj, resolution = 1, use_for_inference = NULL) {
  coords = data.frame(coords)
  scale_factor = zone_obj$scale_factor
  if (scale_factor > 1) {
    # Move into a reasonable range for interpolation
    coords$x = coords$x * scale_factor
    coords$y = coords$y * scale_factor
    resolution = resolution * scale_factor
  }
  if (length(use_for_inference)) {
    coords_subset = coords[rownames(coords) %in% use_for_inference,]
  } else {
    coords_subset = coords
  }
  interp_data = apply_interpolation(coords_subset, zone_obj, resolution)
  ix = findInterval(coords$x, interp_data$x)
  iy = findInterval(coords$y, interp_data$y)
  ix[ix == 0] = 1
  iy[iy == 0] = 1
  ix[ix == length(interp_data$x)] = length(interp_data$x) - 1
  iy[iy == length(interp_data$y)] = length(interp_data$y) - 1
  interp_value = mapply(function(i, j, x0, y0) {
    x1 = interp_data$x[i]; x2 = interp_data$x[i + 1]
    y1 = interp_data$y[j]; y2 = interp_data$y[j + 1]
    z11 = interp_data$z[i, j]; z21 = interp_data$z[i + 1, j]
    z12 = interp_data$z[i, j + 1]; z22 = interp_data$z[i + 1, j + 1]
    if (any(is.na(c(z11, z21, z12, z22)))) return(NA)
    wx = (x0 - x1) / (x2 - x1)
    wy = (y0 - y1) / (y2 - y1)
    z = (1 - wx) * (1 - wy) * z11 + wx * (1 - wy) * z21 +
      (1 - wx) * wy * z12 + wx * wy * z22
    return(z)
  }, ix, iy, coords$x, coords$y)
  coords$zonation = interp_value
  breaks = seq(1, 3, length.out = 4)
  zone = cut(coords$zonation, breaks = breaks, labels = c('Zone_1', 'Zone_2', 'Zone_3'), include.lowest = TRUE)
  zone = factor(zone, levels = ZONES)
  names(zone) = rownames(coords)
  zone
}

#' Plots the 2-dimensional interpolated zones in a spatial dataset
#'
#' @param coords Coordinate matrix with samples as rows, and columns `x` and `y`. Rownames of coords should match colnames of zone_obj$mtx.
#' @param zone_obj Calibrated Zonation Object
#' @param resolution Optional numeric value for the resolution, where higher value results in a more granular interpolation (default 1)
#' @param use_for_inference (optional) A vector of sample names which should be used for zonation inference (recommended to use only hepatocytes, if annotation is available). If not provided, all samples will be used.
#' @return A ggplot object
#' @export
plotZoneSpatial = function(coords, zone_obj, resolution = 1, use_for_inference = NULL) {
  coords = data.frame(coords)
  scale_factor = zone_obj$scale_factor
  if (scale_factor > 1) {
    # Move into a reasonable range for interpolation
    coords$x = coords$x * scale_factor
    coords$y = coords$y * scale_factor
    resolution = resolution * scale_factor
  }
  if (length(use_for_inference)) {
    coords_subset = coords[rownames(coords) %in% use_for_inference,]
  } else {
    coords_subset = coords
  }
  interp_data = apply_interpolation(coords_subset, zone_obj, resolution)
  interp_df = data.frame(
    x = rep(interp_data$x, times = length(interp_data$y)),
    y = rep(interp_data$y, each = length(interp_data$x)),
    z = as.vector(interp_data$z)
  )
  if (scale_factor > 1) {
    # Transfer back to original space
    interp_df$x = interp_df$x / scale_factor
    interp_df$y = interp_df$y / scale_factor
  }
  breaks = seq(1, 3, length.out = 4)
  ggplot(interp_df, aes(x, y, z = z)) +
    geom_contour_filled(breaks = breaks) +
    coord_fixed() +
    scale_fill_manual(name = 'Zone', labels = paste('Zone', 1:3), values = c('#333F48', '#C6AA76', '#BA0C2F')) # "Bull City" color scheme from https://r-graph-gallery.com/color-palette-finder
}

#' Plots the 2-dimensional interpolated zones with zonation contour outlines in a spatial dataset
#'
#' @param coords Coordinate matrix with samples as rows, and columns `x` and `y`. Rownames of coords should match colnames of zone_obj$mtx.
#' @param zone_obj Calibrated Zonation Object
#' @param resolution Optional numeric value for the resolution, where higher value results in a more granular interpolation (default 1)
#' @param point_size Optional numeric value for the ggplot point size (default 1)
#' @param line_width Optional numeric value for the ggplot contour line width (default 2)
#' @param use_for_inference (optional) A vector of sample names which should be used for zonation inference (recommended to use only hepatocytes, if annotation is available). If not provided, all samples will be used.
#' @return A ggplot object
#' @export
plotZoneSpatialContours = function(coords, zone_obj, resolution = 1, point_size = 1, line_width = 2, plot_options = NULL, use_for_inference = NULL) {
  coords = data.frame(coords)
  scale_factor = zone_obj$scale_factor
  if (scale_factor > 1) {
    # Move into a reasonable range for interpolation
    coords$x = coords$x * scale_factor
    coords$y = coords$y * scale_factor
    resolution = resolution * scale_factor
  }
  if (length(use_for_inference)) {
    coords_subset = coords[rownames(coords) %in% use_for_inference,]
  } else {
    coords_subset = coords
  }
  interp_data = apply_interpolation(coords_subset, zone_obj, resolution)
  coords$zone = getZonationGradient(zone_obj)[rownames(coords)]
  interp_df = data.frame(
    x = rep(interp_data$x, times = length(interp_data$y)),
    y = rep(interp_data$y, each = length(interp_data$x)),
    z = as.vector(interp_data$z)
  )
  if (scale_factor > 1) {
    # Transfer back to original space
    coords$x = coords$x / scale_factor
    coords$y = coords$y / scale_factor
    interp_df$x = interp_df$x / scale_factor
    interp_df$y = interp_df$y / scale_factor
  }
  breaks = seq(1, 3, length.out = 4)
  ggplot(coords) +
    geom_point(data = coords, aes(x = x, y = y, color = zone), size = point_size) +
    scale_color_viridis_c(name = 'Zonation') +
    geom_contour(data = interp_df,
                aes(x = x, y = y, z = z),
                breaks = breaks,
                color = 'black',
                linewidth = line_width) +
    coord_fixed() +
    dark_theme_classic()
}

#' Plots a custom variable with zonation contour outlines in a spatial dataset
#'
#' @param meta Metadata matrix with samples as rows, and columns `x`, `y`, and `mycolname`, where `mycolname` is passed as `colname`. Rownames of coords should match colnames of zone_obj$mtx.
#' @param col_name Name of custom column in `meta`
#' @param zone_obj Calibrated Zonation Object
#' @param resolution Optional numeric value for the resolution, where higher value results in a more granular interpolation (default 1)
#' @param point_size Optional numeric value for the ggplot point size (default 1)
#' @param use_for_inference (optional) A vector of sample names which should be used for zonation inference (recommended to use only hepatocytes, if annotation is available). If not provided, all samples will be used.
#' @return A ggplot object
#' @export
plotZoneSpatialCustom = function(meta, col_name, zone_obj, resolution = 1, point_size = 1, use_for_inference = NULL) {
  meta = data.frame(meta)
  scale_factor = zone_obj$scale_factor
  if (scale_factor > 1) {
    # Move into a reasonable range for interpolation
    meta$x = meta$x * scale_factor
    meta$y = meta$y * scale_factor
    resolution = resolution * scale_factor
  }
  if (length(use_for_inference)) {
    meta_subset = meta[rownames(meta) %in% use_for_inference,]
  } else {
    meta_subset = meta
  }
  interp_data = apply_interpolation(meta_subset, zone_obj, resolution)
  interp_df = data.frame(
    x = rep(interp_data$x, times = length(interp_data$y)),
    y = rep(interp_data$y, each = length(interp_data$x)),
    z = as.vector(interp_data$z)
  )
  if (scale_factor > 1) {
    # Transfer back to original space
    meta$x = meta$x / scale_factor
    meta$y = meta$y / scale_factor
    interp_df$x = interp_df$x / scale_factor
    interp_df$y = interp_df$y / scale_factor
  }
  breaks = seq(1, 3, length.out = 4)
  ggplot(meta) +
    geom_point(data = meta, aes(x = x, y = y, color = .data[[col_name]]), size = point_size) +
    geom_contour(data = interp_df,
                aes(x = x, y = y, z = z),
                breaks = breaks,
                color = 'black',
                linewidth = 2) +
    coord_fixed()
}

#' Plot a virtual lobule of the inferred zonation distribution
#'
#' Visualizes the distribution of inferred zonation within an idealized
#' hexagonal lobule. Each pixel is placed according to its rank along the
#' corner-to-center axis (0 = nearest corner / Zone 1, 1 = center / Zone 3)
#' and colored by the matching quantile of \code{getZonationGradient(zone_obj)}.
#'
#' @param zone_obj A calibrated Zonation Object (output of \code{applyModel}).
#' @param resolution Integer pixel resolution along the x-axis (default 100).
#' @param palette Character; a \pkg{paletteer} continuous palette name
#' @param reverse_palette Logical; reverse palette direction (default FALSE).
#' @param pointy_top Logical; if TRUE the hexagon has a vertex at the top,
#'   otherwise it is flat-topped (default FALSE, matching \code{virtual_lobule}).
#' @param show_legend Logical; show the colour legend (default TRUE).
#' @param seed Optional integer for reproducibility. Not strictly required
#'   here (the mapping is deterministic), but kept for API symmetry with
#'   \code{virtual_lobule}.
#' @return A ggplot object.
#' @export
plotVirtualLobule <- function(zone_obj,
                              resolution = 100,
                              palette = "ggthemes::Classic Red-Blue",
                              reverse_palette = T,
                              pointy_top = FALSE,
                              show_legend = TRUE,
                              seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  # --- 0. Pull the inferred zonation gradient (values on [1, 3]) ---
  zonation <- getZonationGradient(zone_obj)
  zonation <- zonation[is.finite(zonation)]
  if (length(zonation) == 0) {
    stop("No finite zonation values from getZonationGradient(zone_obj).")
  }

  # --- 1. Hexagon geometry ---
  R <- 1
  angle_offset <- if (pointy_top) pi / 6 else 0
  vertex_angles <- angle_offset + (0:5) * pi / 3
  vertices_x <- R * cos(vertex_angles)
  vertices_y <- R * sin(vertex_angles)

  # --- 2. Pixel grid ---
  x_range <- range(vertices_x)
  y_range <- range(vertices_y)
  pixel_size <- diff(x_range) / resolution

  x_seq <- seq(x_range[1] + pixel_size / 2, x_range[2] - pixel_size / 2, by = pixel_size)
  y_seq <- seq(y_range[1] + pixel_size / 2, y_range[2] - pixel_size / 2, by = pixel_size)
  grid <- expand.grid(x = x_seq, y = y_seq)

  # --- 3. Point-in-hexagon test (cross-product against CCW edges) ---
  inside <- rep(TRUE, nrow(grid))
  for (i in 1:6) {
    j <- (i %% 6) + 1
    ex <- vertices_x[j] - vertices_x[i]
    ey <- vertices_y[j] - vertices_y[i]
    cross <- ex * (grid$y - vertices_y[i]) - ey * (grid$x - vertices_x[i])
    inside <- inside & (cross >= 0)
  }
  grid <- grid[inside, , drop = FALSE]
  if (nrow(grid) == 0) stop("No pixels inside hexagon. Increase resolution.")

  # --- 4. Position: 0 = nearest corner (portal), 1 = center (central) ---
  d_center <- sqrt(grid$x^2 + grid$y^2)
  d_nearest_corner <- do.call(pmin, lapply(seq_len(6), function(i) {
    sqrt((grid$x - vertices_x[i])^2 + (grid$y - vertices_y[i])^2)
  }))
  raw_pos <- d_nearest_corner / (d_nearest_corner + d_center)
  raw_pos[d_center < 1e-12] <- 1

  # Rank-normalise pixel positions to a uniform [0, 1] distribution
  r <- rank(raw_pos, ties.method = "average")
  grid$position <- (r - min(r)) / (max(r) - min(r))

  # --- 5. Map each pixel to the matching quantile of the zonation gradient ---
  grid$zonation <- as.numeric(stats::quantile(
    zonation,
    probs = grid$position,
    na.rm = TRUE,
    names = FALSE,
    type  = 7
  ))

  # --- 6. Colour palette ---
  n_colors <- 256
  pal_colors <- tryCatch(
    as.character(paletteer::paletteer_c(palette, n = n_colors)),
    error = function(e) as.character(paletteer::paletteer_d(palette))
  )
  if (reverse_palette) pal_colors <- rev(pal_colors)

  # --- 7. ggplot ---
  ggplot2::ggplot(grid, ggplot2::aes(x = x, y = y, fill = zonation)) +
    ggplot2::geom_tile(width = pixel_size, height = pixel_size) +
    ggplot2::scale_fill_gradientn(
      colours  = pal_colors,
      limits   = c(1, 3),
      oob      = scales::squish,
      na.value = "grey50",
      name     = "Zonation"
    ) +
    ggplot2::coord_fixed(expand = FALSE) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position  = if (show_legend) "right" else "none",
      plot.background  = ggplot2::element_rect(fill = "transparent", colour = NA),
      panel.background = ggplot2::element_rect(fill = "transparent", colour = NA)
    )
}

#' Density plot of inferred zonation across one or more samples
#'
#' Overlays the distributions of \code{getZonationGradient()} values for one
#' or more samples as outlined density curves on a shared baseline.
#'
#' @param zone_objs A single Zonation Object, or a (preferably named) list of
#'   them. List names are used as legend labels.
#' @param palette Character; a paletteer palette name for outline colors
#'   (default \code{"grDevices::rainbow"}).
#' @param line_width Numeric line width for density outlines (default 1).
#' @param adjust Bandwidth adjustment passed to \code{geom_density} (default 3).
#' @param log_y Logical; if TRUE, plot density on a log y-scale (default FALSE).
#' @return A ggplot object.
#' @export
plotZonationRidge <- function(zone_objs, palette = 'grDevices::rainbow', line_width = 1, adjust = 3, log_y = FALSE) {
  if (is.list(zone_objs) && !is.null(zone_objs$mtx)) {
    zone_objs <- list(zone_objs)
  }
  if (is.null(names(zone_objs))) {
    names(zone_objs) <- paste0("Sample_", seq_along(zone_objs))
  } else {
    blank <- !nzchar(names(zone_objs))
    names(zone_objs)[blank] <- paste0("Sample_", which(blank))
  }
  df <- do.call(rbind, lapply(seq_along(zone_objs), function(i) {
    z <- as.numeric(getZonationGradient(zone_objs[[i]]))
    z <- z[is.finite(z)]
    if (length(z) == 0) {
      stop(sprintf("No finite zonation values for sample '%s'.",
                   names(zone_objs)[i]))
    }
    data.frame(Sample = names(zone_objs)[i], Zonation = z)
  }))
  df$Sample <- factor(df$Sample, levels = names(zone_objs))
  n_samples <- length(zone_objs)
  pal_colors <- tryCatch({
    cols <- as.character(paletteer::paletteer_d(palette))
    rep_len(cols, n_samples)
  }, error = function(e) {
    as.character(paletteer::paletteer_c(palette, n = n_samples))
  })

  # Build the base plot
  bw_fixed <- stats::bw.nrd0(df$Zonation)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = Zonation, color = Sample)) +
    ggplot2::geom_density(
      fill = NA,
      linewidth = line_width,
      bounds = c(1, 3),
      adjust = adjust,
      n = 1024,
      bw = bw_fixed
    ) +
    ggplot2::scale_color_manual(values = pal_colors) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::theme_gray()

  # Conditionally apply log scale or standard linear scale
  if (log_y) {
    p <- p +
      ggplot2::scale_y_log10() +
      # The ylim here acts as a floor to hide the microscopic floating-point noise
      ggplot2::coord_cartesian(xlim = c(1, 3), ylim = c(1e-4, NA)) +
      ggplot2::labs(x = "Zonation", y = "Log Proportion")
  } else {
    p <- p +
      ggplot2::coord_cartesian(xlim = c(1, 3)) +
      ggplot2::labs(x = "Zonation", y = "Proportion")
  }

  return(p)
}

#' Compare per-gene zonation between two samples
#'
#' Heatmap of each landmark gene's posterior layer distribution from the
#' baseline (sample x) fit, with a side bar showing the log fold-change of
#' mean fraction from x to y. Rows are the top zone-1 and zone-3 reference
#' genes ranked by how concentrated their posterior layer distribution is
#' toward the periportal or pericentral end of the lobule.
#'
#' @param zone_obj_x Reference (baseline) zonation object.
#' @param zone_obj_y Comparison zonation object (must share the same fit as x).
#' @param n_per_zone Number of top genes per zone to display. Default 30.
#' @param font_size Axis label font size. Default 9.
#' @return A patchwork object (heatmap + side bar).
#' @export
plotZonationHeat = function(zone_obj_x, zone_obj_y, n_per_zone = 30, font_size = 9) {
  if (!requireNamespace('patchwork', quietly = TRUE))
    stop("patchwork is required; install with install.packages('patchwork').")
  if (!identical(zone_obj_x$fit$lm, zone_obj_y$fit$lm) ||
      zone_obj_x$fit$r0 != zone_obj_y$fit$r0)
    stop('zone_obj_x and zone_obj_y must share the same fit.')

  mtx_x = zone_obj_x$mtx
  mtx_y = zone_obj_y$mtx
  n_layers = zone_obj_x$fit$n_layers
  eps = 1e-6

  m_post = zone_obj_x$fit$alpha / zone_obj_x$fit$beta
  row_tot = rowSums(m_post, na.rm = TRUE)
  row_tot[!is.finite(row_tot) | row_tot < 1e-12] = 1e-12
  p_layer = sweep(m_post, 1, row_tot, '/')

  common = intersect(rownames(p_layer), rownames(mtx_y))
  if (length(common) == 0) stop('no reference genes present in y.')
  p_layer = p_layer[common, , drop = FALSE]

  Nx = Matrix::colSums(mtx_x); Nx[Nx == 0] = 1
  Ny = Matrix::colSums(mtx_y); Ny[Ny == 0] = 1
  mu_x = rowMeans(sweep(as.matrix(mtx_x[common, , drop = FALSE]), 2, Nx, '/'))
  mu_y = rowMeans(sweep(as.matrix(mtx_y[common, , drop = FALSE]), 2, Ny, '/'))
  log_expr = log1p(mu_x * 1e4)

  midpoint = (n_layers + 1) / 2
  com = as.numeric(p_layer %*% seq_len(n_layers))
  names(com) = rownames(p_layer)
  zonation = (com - midpoint) * log_expr

  n_z1 = n_z3 = n_per_zone
  z1_pool = names(sort(zonation[is.finite(zonation) & zonation < 0]))
  z3_pool = names(sort(zonation[is.finite(zonation) & zonation > 0], decreasing = TRUE))
  z1_top = head(z1_pool, n_z1)
  z3_top = head(z3_pool, n_z3)
  ordered_genes = c(z1_top, rev(z3_top))
  if (length(ordered_genes) == 0) stop('no genes after zone filtering.')

  log_fc_top = log((mu_y[ordered_genes] + eps) / (mu_x[ordered_genes] + eps))

  p_layer_top = p_layer[ordered_genes, , drop = FALSE]

  heat_df = data.frame(
    gene  = factor(rep(ordered_genes, times = n_layers), levels = rev(ordered_genes)),
    layer = factor(rep(seq_len(n_layers), each = length(ordered_genes)),
                   levels = seq_len(n_layers)),
    value = as.numeric(p_layer_top)
  )
  bar_df = data.frame(
    gene   = factor(ordered_genes, levels = rev(ordered_genes)),
    log_fc = log_fc_top
  )

  fill_max = max(heat_df$value, na.rm = TRUE); if (!is.finite(fill_max) || fill_max == 0) fill_max = 1
  bar_lim  = max(abs(bar_df$log_fc), na.rm = TRUE); if (!is.finite(bar_lim) || bar_lim == 0) bar_lim = 1

  p_heat = ggplot2::ggplot(heat_df, ggplot2::aes(x = layer, y = gene, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(limits = c(0, fill_max), na.value = 'grey90',
                                  name = 'Model\nPosterior') +
    ggplot2::labs(x = 'Layer', y = 'Gene') +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text = ggplot2::element_text(size = font_size))

  p_bar = ggplot2::ggplot(bar_df, ggplot2::aes(x = log_fc, y = gene, fill = log_fc > 0)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = c('TRUE' = '#b40426', 'FALSE' = '#3b4cc0')) +
    ggplot2::scale_x_continuous(limits = c(-bar_lim, bar_lim)) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::labs(x = 'log FC\n(y vs x)') +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                   axis.text.y  = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   legend.position = 'none',
                   axis.text = ggplot2::element_text(size = font_size))

  patchwork::wrap_plots(p_heat, p_bar, widths = c(4, 1))
}
