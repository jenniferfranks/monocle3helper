#' Plot a spatial signature overlay for a single Xenium sample
#'
#' Overlays a numeric signature score stored in colData(cds) onto a
#' Xenium morphology image, with cell centroids converted into the
#' same pixel space as the image returned by \code{build_xenium_images()}.
#'
#' @param cds A monocle3 cell_data_set containing cell centroid
#'   coordinates and the signature column.
#' @param images List returned by build_xenium_images().
#' @param signature_col Name of numeric column in colData(cds).
#' @param sample_id Sample identifier matching colData(cds)$sample.
#' @param point_size Size of cell centroid points.
#' @param point_alpha Transparency of cell centroid points.
#' @param img_buffer_only_cells Logical; crop the view to the region
#'   containing cells rather than showing the full image extent.
#' @param palette Color palette type.
#' @param custom_cols Custom colors if palette = "custom".
#' @param viridis_opt Viridis palette option.
#' @param diverging_midpoint Midpoint for diverging color scale.
#'
#' @return A ggplot object.
#'
#' @export
plot_signature_spatial_sample_xenium <- function(
    cds,
    images,
    signature_col,
    sample_id,
    point_size  = 0.3,
    point_alpha = 0.8,
    img_buffer_only_cells = TRUE,
    palette      = c("diverging", "sequential", "viridis", "custom"),
    custom_cols  = NULL,
    viridis_opt  = "magma",
    diverging_midpoint = 0
) {

  palette <- match.arg(palette)

  ## -------------------- sanity checks --------------------
  cd <- SummarizedExperiment::colData(cds)

  if (!signature_col %in% colnames(cd)) {
    stop("signature_col not found in colData(cds): ", signature_col)
  }

  if (!is.numeric(cd[[signature_col]])) {
    stop("signature_col must be numeric: ", signature_col)
  }

  df_cells <- as.data.frame(cd)
  df_cells <- df_cells[df_cells$sample == sample_id, , drop = FALSE]

  if (nrow(df_cells) == 0) {
    stop("No cells found for sample_id: ", sample_id)
  }

  required <- c("x_centroid", "y_centroid")
  if (!all(required %in% colnames(df_cells))) {
    stop("Cell centroid coordinates missing from colData(cds).")
  }

  ## -------------------- image --------------------
  img_entry <- images$images[[sample_id]]
  if (is.null(img_entry)) {
    stop(
      "No morphology image available for sample_id: ", sample_id,
      ". Did you run build_xenium_images()?"
    )
  }

  ## -------------------- centroid coordinates in image pixel space --------------------
  effective_pixel_size <- images$effective_pixel_size_um
  df_cells$x_pixel <- df_cells$x_centroid / effective_pixel_size
  df_cells$y_pixel <- df_cells$y_centroid / effective_pixel_size

  view_xlim <- if (img_buffer_only_cells) img_entry$cells_xlim else img_entry$xlim
  view_ylim <- if (img_buffer_only_cells) img_entry$cells_ylim else img_entry$ylim

  xy_ratio <- diff(view_xlim) / diff(view_ylim)

  ## -------------------- color scale --------------------
  color_scale <- switch(
    palette,
    "sequential" = ggplot2::scale_colour_gradient(
      low = "white", high = "red", na.value = "grey80"
    ),
    "diverging" = ggplot2::scale_colour_gradient2(
      low = "blue",
      mid = "grey90",
      high = "red",
      midpoint = diverging_midpoint,
      na.value = "grey80"
    ),
    "viridis" = ggplot2::scale_colour_viridis_c(
      option = viridis_opt,
      na.value = "grey80"
    ),
    "custom" = {
      if (is.null(custom_cols) || length(custom_cols) < 2) {
        stop("custom_cols must contain at least two colors.")
      }
      ggplot2::scale_colour_gradientn(
        colors = custom_cols,
        na.value = "grey80"
      )
    }
  )

  sample_label <- unique(df_cells$sample_label)

  ## -------------------- plot --------------------
  # annotation_raster() always draws row 1 of the raster at whatever data-y
  # value is passed as ymax; ymin/ymax are swapped here (relative to the
  # image's natural 0..height extent) so that, combined with
  # scale_y_reverse() below, the image keeps the same orientation as cell
  # points (whose y_pixel also increases downward, matching image rows).
  ggplot2::ggplot() +
    ggplot2::annotation_raster(
      img_entry$raster,
      xmin = img_entry$xlim[1], xmax = img_entry$xlim[2],
      ymin = img_entry$ylim[2], ymax = img_entry$ylim[1]
    ) +
    ggplot2::coord_cartesian(xlim = view_xlim, ylim = view_ylim, expand = FALSE) +
    ggplot2::scale_y_reverse() +
    ggplot2::geom_point(
      data = df_cells,
      ggplot2::aes(
        x = x_pixel,
        y = y_pixel,
        colour = .data[[signature_col]]
      ),
      size  = point_size,
      alpha = point_alpha
    ) +
    color_scale +
    ggplot2::theme_void() +
    ggplot2::theme(
      aspect.ratio = 1 / xy_ratio,
      legend.text  = ggplot2::element_text(size = 10, color = "black"),
      legend.title = ggplot2::element_text(size = 10, color = "black")
    ) +
    ggplot2::labs(
      title  = paste0("Sample: ", sample_label, " (", sample_id, ")"),
      colour = signature_col
    )
}
