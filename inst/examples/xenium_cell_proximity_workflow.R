# Example: Xenium cell-type spatial proximity analysis
#
# This assumes `cds` is a Xenium monocle3 cell_data_set created with
# build_xenium_cds() and that cell type labels have been added to colData(cds)
# in a column named "cell_type".

library(monocle3helper)

proximity <- analyze_xenium_cell_proximity(
  cds,
  cell_type_col = "cell_type",
  strata_cols = "sample",
  k = 10,
  max_distance = 50
)

# Directed focal-type -> neighbor-type summaries with abundance-adjusted
# enrichment ratios and distance summaries.
proximity$pairwise

# Cell type abundance and retained-neighbor counts per sample or ROI stratum.
proximity$by_stratum

# One row per retained focal-cell -> neighbor-cell relationship.
proximity$cell_neighbors

plot_xenium_cell_proximity(proximity, value = "enrichment")
plot_xenium_cell_proximity(proximity, value = "median_distance")

# To analyze ROIs independently, include an ROI column created by
# subset_xenium_rois() or another ROI assignment workflow.
roi_proximity <- analyze_xenium_cell_proximity(
  cds,
  cell_type_col = "cell_type",
  strata_cols = "sample",
  roi_col = "roi_id",
  k = 10,
  max_distance = 50,
  exclude_same_type = TRUE
)

# Compare directed proximity statistics between samples. The pairwise table
# keeps the sample column because sample was included in strata_cols.
sample_differences <- compare_xenium_cell_proximity(
  proximity,
  stratum_col = "sample",
  reference = "control_sample",
  value = "enrichment"
)

# If both sample and ROI are in strata_cols, compare samples within matched ROIs.
roi_sample_differences <- compare_xenium_cell_proximity(
  roi_proximity,
  stratum_col = "sample",
  within_cols = "roi_id",
  reference = "control_sample",
  value = "median_distance"
)
