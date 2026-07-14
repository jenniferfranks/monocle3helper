# Example: build and test a Xenium ROI-only CDS
#
# This script shows two ways to test the ROI workflow:
#   1. build a smaller ROI-only CDS directly from Xenium output files; or
#   2. subset an already-built Xenium CDS to the same ROI definitions.
#
# Edit the "User inputs" section before running.

## -------------------------------------------------------------------------
## Setup
## -------------------------------------------------------------------------

# If you are running from a checkout of monocle3helper, load the local package.
# Otherwise, replace this with library(monocle3helper).
if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
  devtools::load_all(".")
} else {
  library(monocle3helper)
}

library(SummarizedExperiment)

## -------------------------------------------------------------------------
## User inputs: edit these values for your Xenium run
## -------------------------------------------------------------------------

# Directory containing one subdirectory per Xenium sample/run.
xenium_root <- "/path/to/xenium_outputs"

# Optional: specify sample directories and labels explicitly. If sample_dirs is
# NULL, discover_xenium_sample_table() will use all first-level directories
# under xenium_root.
sample_dirs <- NULL
sample_ids <- NULL
sample_labels <- NULL
sample_groups <- NULL

# Output files for quick inspection/reuse.
out_dir <- "xenium_roi_test_outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## -------------------------------------------------------------------------
## Discover Xenium files
## -------------------------------------------------------------------------

sample_table <- discover_xenium_sample_table(
  root_dir = xenium_root,
  sample_dirs = sample_dirs,
  sample_ids = sample_ids,
  sample_labels = sample_labels,
  sample_groups = sample_groups
)

print(sample_table)

## -------------------------------------------------------------------------
## Automatically derive smaller ROIs per sample
## -------------------------------------------------------------------------

# This reads each cells.csv.gz listed in sample_table and computes a rectangular
# ROI around the observed Xenium cell centroids, using the same centroid-range
# idea that build_xenium_images() uses for its cell-focused image crop.
# This example keeps about 10% of each sample extent across 5 reproducible,
# randomly placed, non-overlapping ROIs.
roi_table <- derive_xenium_roi_table(
  sample_table = sample_table,
  n_rois_per_sample = 5,
  target_fraction = 0.1,
  selection = "random",
  random_seed = 12345,
  allow_overlap = FALSE,
  roi_id_prefix = "auto_roi"
)

print(roi_table)
write.csv(
  roi_table,
  file.path(out_dir, "auto_roi_table.csv"),
  row.names = FALSE
)

## -------------------------------------------------------------------------
## Option A: build a small ROI-only CDS directly
## -------------------------------------------------------------------------

cds_roi <- build_xenium_cds(
  sample_table = sample_table,
  roi_table = roi_table,
  roi_id_col = "roi_id",
  roi_sample_col = "sample",
  x_min_col = "x_min",
  x_max_col = "x_max",
  y_min_col = "y_min",
  y_max_col = "y_max",
  num_dim = 30,
  k = 15,
  random_seed = 12345
)

saveRDS(cds_roi, file.path(out_dir, "xenium_roi_only_cds.rds"))

roi_summary <- as.data.frame(colData(cds_roi)) |>
  dplyr::count(sample, roi_id, name = "n_cells")

print(roi_summary)
write.csv(
  roi_summary,
  file.path(out_dir, "roi_cell_counts.csv"),
  row.names = FALSE
)

## -------------------------------------------------------------------------
## Option B: subset an existing/full CDS to the same ROIs
## -------------------------------------------------------------------------

# If you already have a full Xenium CDS, load it here and run this block.
# full_cds <- readRDS("/path/to/full_xenium_cds.rds")
# cds_roi_from_full <- subset_xenium_rois(
#   cds = full_cds,
#   roi_table = roi_table,
#   roi_id_col = "roi_id",
#   roi_sample_col = "sample",
#   cds_sample_col = "sample"
# )
# saveRDS(cds_roi_from_full, file.path(out_dir, "xenium_roi_from_full_cds.rds"))

## -------------------------------------------------------------------------
## Optional quick spatial check with ROI cell centroids
## -------------------------------------------------------------------------

roi_cells <- as.data.frame(colData(cds_roi))

p <- ggplot2::ggplot(
  roi_cells,
  ggplot2::aes(x = x_centroid, y = y_centroid, color = roi_id)
) +
  ggplot2::geom_point(size = 0.2, alpha = 0.7) +
  ggplot2::coord_fixed() +
  ggplot2::facet_wrap(~ sample) +
  ggplot2::theme_minimal() +
  ggplot2::labs(
    title = "Xenium ROI-only CDS cell centroids",
    x = "x centroid (microns)",
    y = "y centroid (microns)",
    color = "ROI"
  )

print(p)
ggplot2::ggsave(
  filename = file.path(out_dir, "roi_centroid_check.png"),
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)
