plot_umap_manustyle <- function(
  cds,
  celltype_colors = NULL,
  color_cells_by = "fine_annotations",
  group_cells_by = "fine_annotations",
  palette = scales::hue_pal(),
  cell_size = 0.1,
  label_size = 1.5,
  label_fill_alpha = 0.7,
  show_trajectory_graph = FALSE,
  label_cell_groups = FALSE
) {

  # Extract labels
  labels <- as.character(colData(cds)[[group_cells_by]])

  # Generate colors if not supplied
  if (is.null(celltype_colors)) {
    celltypes <- sort(unique(labels))
    celltype_colors <- setNames(
      palette(length(celltypes)),
      celltypes
    )
  }

  # Base Monocle plot
  p <- plot_cells(
    cds,
    color_cells_by = color_cells_by,
    group_cells_by = group_cells_by,
    label_cell_groups = label_cell_groups,
    cell_size = cell_size,
    show_trajectory_graph = show_trajectory_graph
  )

  # Extract UMAP coordinates
  df <- as.data.frame(reducedDims(cds)$UMAP)
  colnames(df) <- c("x", "y")
  df$label <- labels

  # Compute centroids
  centroids <- df |>
    dplyr::group_by(label) |>
    dplyr::summarize(
      x = mean(x),
      y = mean(y),
      .groups = "drop"
    )

  p +
    ggplot2::scale_color_manual(
      values = celltype_colors,
      na.value = "#F5F5F5"
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none") +
    ggrepel::geom_label_repel(
      data = centroids,
      ggplot2::aes(x, y, label = label, color = label),
      size = label_size,
      box.padding = 0.2,
      point.padding = 0.3,
      label.padding = grid::unit(0.1, "lines"),
      label.r = grid::unit(0.1, "lines"),
      fill = scales::alpha("white", label_fill_alpha),
      segment.color = "grey70",
      min.segment.length = 0
    )
}