Package dedicated to storing common helper functions for Franks lab


remotes::install_github("jenniferfranks/monocle3helper")

## Xenium ROI downsampling and POLC proximity workflow

For very large Xenium case/control samples, start with reproducible spatial
regions of interest (ROIs) rather than random cells. This preserves local tissue
neighborhoods and makes the downsample traceable.

```r
# Pick 5-10 ROIs per sample. Tune roi_width_um/roi_height_um so each ROI is
# large enough to contain POLCs plus neighbors but small enough for monocle3.
rois <- select_xenium_rois(
  cds,
  n_rois = 8,
  roi_width_um = 1000,
  roi_height_um = 1000,
  strategy = "quantile_grid",
  min_cells = 500,
  seed = 20260714
)

# Save this table with the analysis outputs so every downsample is auditable.
write.csv(rois, "xenium_roi_table.csv", row.names = FALSE)

# Plot where each ROI came from in the original sample coordinate system.
p <- plot_xenium_rois(cds, rois)
ggplot2::ggsave("xenium_roi_overview.png", p, width = 10, height = 6, dpi = 300)

# Subset to cells inside the selected ROIs, then run the heavier monocle3 steps.
cds_roi <- subset_xenium_rois(cds, rois)

# With broad cell type annotations in colData(cds_roi)$broad_cell_type, ask which
# cell type is the nearest neighbor of each POLC. Set exclude_same_cell_type = TRUE
# if the biological question is specifically about non-POLC neighbors.
polc_nn <- nearest_cell_type_summary(
  cds_roi,
  focal_types = "POLC",
  cell_type_col = "broad_cell_type",
  exclude_same_cell_type = TRUE
)
polc_nn$summary
```

Recommended strategy:

1. Choose the same ROI size and number of ROIs for the case and control sample.
2. Use `strategy = "quantile_grid"` first so ROIs span low-, medium-, and
   high-cell-density tissue regions instead of only the densest area.
3. Save the ROI table and overview figure with every run; the ROI bounds are in
   Xenium micron coordinates and can be reused exactly.
4. Run POLC proximity summaries within the ROI subset and compare neighbor-cell
   proportions by sample, condition, and ROI. For robustness, repeat the ROI
   selection with a different seed or increase `n_rois` if findings are driven
   by a single ROI.
