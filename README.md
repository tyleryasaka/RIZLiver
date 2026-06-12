# RIZLiver
Robust inference of liver zonation from gene expression

If you have any questions or comments, or any issues using the package, please feel free to reach out: `typublic@pm.me`

## Installation

```r
install.packages("devtools")
devtools::install_github("tyleryasaka/RIZLiver")
```

## Usage

### Obtaining zonation in Seurat v5

```r
library(RIZLiver)

# Train the model on a baseline sample
seurat_baseline <- subset(seurat_obj, subset = sample == 'my_baseline_sample_id')
zonation_obj <- trainModel(Seurat::GetAssayData(seurat_baseline, layer = 'counts'), species = 'human')

# Apply the trained model to your full dataset (or any other sample)
zonation_obj <- applyModel(Seurat::GetAssayData(seurat_obj, layer = 'counts'), zonation_obj)

# Obtain discrete zonation bins, normalized to the baseline
zonation_assignments <- getZone(zonation_obj)
seurat_obj <- AddMetaData(seurat_obj, zonation_assignments, col.name = 'zone')

# Obtain a continuous zonation gradient (values between 1 and 3), normalized to the baseline
zonation_gradient <- getZonationGradient(zonation_obj)
seurat_obj <- AddMetaData(seurat_obj, zonation_gradient, col.name = 'zonation')
```

## Functions

### Model setup

#### `trainModel()`

Train the model on a baseline liver sample.

**Usage:**
```r
trainModel(mtx, coords = NULL, species = 'human', verbose = FALSE)
```

**Arguments:**
- `mtx`: Gene expression matrix (*raw counts*) with genes as rows.
- `coords`: (Optional) Coordinate matrix with samples as rows and columns `x` and `y`. Used to set the spatial scale factor for downstream spatial functions.
- `species`: (Optional) Species to use; defaults to `'human'`. Supports `'human'` and `'mouse'`.
- `verbose`: (Optional) If `TRUE`, prints training diagnostics (default `FALSE`).

**Returns:**
- A `ZonationObject` with calibrated baseline zonation.

---

#### `applyModel()`

Attach a (new) gene expression matrix to a trained `ZonationObject`. Almost every other function in the package operates on the resulting object.

**Usage:**
```r
applyModel(mtx, zone_obj, coords = NULL)
```

**Arguments:**
- `mtx`: Gene expression matrix (*raw counts*) with genes as rows.
- `zone_obj`: Trained `ZonationObject` (output of `trainModel()`).
- `coords`: (Optional) Coordinate matrix with samples as rows and columns `x` and `y`. Used to set the spatial scale factor for downstream interpolation.

**Returns:**
- A `ZonationObject` with `mtx` attached, ready for downstream queries and plotting.

---

### Extracting zonation

#### `getZone()`

Return the discrete zonation bin for each cell/spot.

**Usage:**
```r
getZone(zone_obj)
```

**Arguments:**
- `zone_obj`: Calibrated `ZonationObject` with a matrix applied (see `applyModel()`).

**Returns:**
- Factor of zonation assignments (`Zone_1`, `Zone_2`, `Zone_3`), one per cell/spot.

---

#### `getZonationGradient()`

Return the continuous zonation value (between 1 and 3) for each cell/spot.

**Usage:**
```r
getZonationGradient(zone_obj)
```

**Arguments:**
- `zone_obj`: Calibrated `ZonationObject` with a matrix applied.

**Returns:**
- Numeric vector of zonation values on `[1, 3]`, one per cell/spot.

---

#### `getGeneZonation()`

Return the per-layer mean expression fraction for each landmark gene, as learned from the reference and updated by the baseline.

**Usage:**
```r
getGeneZonation(zone_obj)
```

**Arguments:**
- `zone_obj`: Calibrated `ZonationObject`.

**Returns:**
- A data frame with one row per landmark gene and one column per layer (`layer_1`, `layer_2`, ...). Values are the mean expression fractions used in the model.

---

### Plotting zonation

#### `plotRegression()`

Scatter and cubic polynomial fit of a gene's expression along the zonation axis.

**Usage:**
```r
plotRegression(zone_obj, gene)
```

**Arguments:**
- `zone_obj`: Calibrated `ZonationObject` with a matrix applied.
- `gene`: Name of the gene to plot (must be a row in the applied matrix).

**Returns:**
- A `ggplot` object.

---

#### `plotVirtualLobule()`

Plot an idealized hexagonal lobule, with each pixel placed by its rank along the corner-to-center axis (0 = nearest corner / Zone 1, 1 = center / Zone 3) and colored by the matching quantile of the inferred zonation gradient.

**Usage:**
```r
plotVirtualLobule(zone_obj,
                  resolution = 100,
                  palette = "ggthemes::Classic Red-Blue",
                  reverse_palette = TRUE,
                  pointy_top = FALSE,
                  show_legend = TRUE,
                  seed = NULL)
```

**Arguments:**
- `zone_obj`: Calibrated `ZonationObject` with a matrix applied.
- `resolution`: (Optional) Integer pixel resolution along the x-axis (default 100).
- `palette`: (Optional) A `paletteer` continuous palette name (default `"ggthemes::Classic Red-Blue"`).
- `reverse_palette`: (Optional) If `TRUE`, reverse the palette direction (default `TRUE`).
- `pointy_top`: (Optional) If `TRUE`, the hexagon has a vertex at the top; otherwise flat-topped (default `FALSE`).
- `show_legend`: (Optional) If `TRUE`, show the color legend (default `TRUE`).
- `seed`: (Optional) Integer seed for reproducibility.

**Returns:**
- A `ggplot` object.

---

#### `plotZonationRidge()`

Overlay the inferred zonation distributions of one or more samples as outlined density curves.

**Usage:**
```r
plotZonationRidge(zone_objs, palette = 'grDevices::rainbow', line_width = 1, adjust = 3, log_y = FALSE)
```

**Arguments:**
- `zone_objs`: A single `ZonationObject`, or a (preferably named) list of them. List names are used as legend labels.
- `palette`: (Optional) A `paletteer` palette name for the outline colors (default `'grDevices::rainbow'`).
- `line_width`: (Optional) Line width for density outlines (default 1).
- `adjust`: (Optional) Bandwidth adjustment passed to `geom_density` (default 3).
- `log_y`: (Optional) If `TRUE`, plot density on a log y-scale (default `FALSE`).

**Returns:**
- A `ggplot` object.

---

### Zonation heatmap

#### `plotZonationHeat()`

Heatmap of each landmark gene's posterior layer distribution from the baseline (sample x) fit, with a side bar showing the log fold-change of mean fraction from x to y. Rows are the top zone-1 and zone-3 reference genes ranked by how concentrated their posterior layer distribution is toward the periportal or pericentral end of the lobule.

**Usage:**
```r
plotZonationHeat(zone_obj_x, zone_obj_y, n_per_zone = 30, font_size = 9)
```

**Arguments:**
- `zone_obj_x`: Reference (baseline) `ZonationObject`.
- `zone_obj_y`: Comparison `ZonationObject` (must share the same fit as `zone_obj_x`).
- `n_per_zone`: (Optional) Number of top genes per zone to display (default 30).
- `font_size`: (Optional) Axis label font size (default 9).

**Returns:**
- A patchwork object combining the heatmap and side bar.

---

### Spatial zonation

These functions infer zonation from gene expression and then propagate it across a tissue section via spatial interpolation. This is particularly useful when you want to assign zonation to non-parenchymal cells based on their proximity to zonated hepatocytes rather than their own expression.

#### `getZoneSpatial()`

Infer discrete zonation for all cells/spots based on spatial interpolation.

**Usage:**
```r
getZoneSpatial(coords, zone_obj, resolution = 1, use_for_inference = NULL)
```

**Arguments:**
- `coords`: Coordinate matrix with samples as rows, and columns `x` and `y`. Rownames should match the colnames of the matrix applied to `zone_obj`.
- `zone_obj`: Calibrated `ZonationObject` with a matrix applied.
- `resolution`: (Optional) Interpolation granularity (default 1; higher = finer).
- `use_for_inference`: (Optional) A vector of sample names which should be used for zonation inference (recommended: hepatocytes only, if annotation is available). If not provided, all samples are used.

**Returns:**
- Factor of discrete zonation assignments for all samples.

---

#### `plotZoneSpatial()`

Plot the 2D interpolated zones as a filled contour map.

**Usage:**
```r
plotZoneSpatial(coords, zone_obj, resolution = 1, use_for_inference = NULL)
```

**Arguments:**
- `coords`: Coordinate matrix with columns `x` and `y`.
- `zone_obj`: Calibrated `ZonationObject` with a matrix applied.
- `resolution`: (Optional) Interpolation granularity (default 1).
- `use_for_inference`: (Optional) Subset of samples for inference.

**Returns:**
- A `ggplot` object with filled contour zones.

---

#### `plotZoneSpatialContours()`

Plot 2D zonation with contour outlines overlaid on the original points.

**Usage:**
```r
plotZoneSpatialContours(coords, zone_obj, resolution = 1, point_size = 1, line_width = 2, plot_options = NULL, use_for_inference = NULL)
```

**Arguments:**
- `coords`: Coordinate matrix with columns `x` and `y`.
- `zone_obj`: Calibrated `ZonationObject` with a matrix applied.
- `resolution`: (Optional) Interpolation granularity (default 1).
- `point_size`: (Optional) Point size (default 1).
- `line_width`: (Optional) Contour line width (default 2).
- `plot_options`: (Optional) Reserved for additional plot customization.
- `use_for_inference`: (Optional) Subset of samples for inference.

**Returns:**
- A `ggplot` object with points colored by zonation and contour outlines.

---

#### `plotZoneSpatialCustom()`

Plot a custom variable spatially, with zonation contour outlines overlaid.

**Usage:**
```r
plotZoneSpatialCustom(meta, col_name, zone_obj, resolution = 1, point_size = 1, use_for_inference = NULL)
```

**Arguments:**
- `meta`: Metadata data frame with samples as rows, and columns `x`, `y`, and the column named by `col_name`. Rownames should match the colnames of the matrix applied to `zone_obj`.
- `col_name`: Name of the column in `meta` to color points by.
- `zone_obj`: Calibrated `ZonationObject` with a matrix applied.
- `resolution`: (Optional) Interpolation granularity (default 1).
- `point_size`: (Optional) Point size (default 1).
- `use_for_inference`: (Optional) Subset of samples for inference.

**Returns:**
- A `ggplot` object with points colored by `col_name` and contour outlines.
