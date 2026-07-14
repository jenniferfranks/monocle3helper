#' Select rectangular regions of interest from large Xenium samples
#'
#' Creates reproducible, per-sample rectangular regions of interest (ROIs) from
#' Xenium cell centroid coordinates. The returned ROI table can be reviewed,
#' plotted with \code{plot_xenium_rois()}, and used to subset a large
#' \code{cell_data_set} with \code{subset_xenium_rois()} before running slower
#' monocle3 steps.
#'
#' @param cds A monocle3 cell_data_set with Xenium centroid coordinates in
#'   \code{colData(cds)}.
#' @param sample_ids Optional character vector of samples to include. Defaults
#'   to all samples in \code{sample_col}.
#' @param n_rois Number of ROIs to choose per sample.
#' @param roi_width_um,roi_height_um ROI width and height in microns.
#' @param strategy ROI placement strategy. \code{"quantile_grid"} scores a grid
#'   of candidate windows and chooses windows spread across cell-density
#'   quantiles; \code{"random"} samples candidate windows at random.
#' @param min_cells Minimum cells required in an ROI candidate.
#' @param sample_col,x_col,y_col Column names in \code{colData(cds)}.
#' @param seed Random seed used for tie-breaking and random sampling.
#' @param grid_n Number of candidate window centers per axis for
#'   \code{"quantile_grid"}.
#'
#' @return A data.frame with one row per selected ROI, including sample, bounds,
#'   center coordinates, and number of cells.
#' @export
select_xenium_rois <- function(
    cds,
    sample_ids = NULL,
    n_rois = 5,
    roi_width_um = 1000,
    roi_height_um = 1000,
    strategy = c("quantile_grid", "random"),
    min_cells = 100,
    sample_col = "sample",
    x_col = "x_centroid",
    y_col = "y_centroid",
    seed = 12345,
    grid_n = 25
) {
  strategy <- match.arg(strategy)
  cd <- as.data.frame(SummarizedExperiment::colData(cds))
  required <- c(sample_col, x_col, y_col)
  if (!all(required %in% colnames(cd))) {
    stop("cds colData must contain: ", paste(required, collapse = ", "))
  }
  if (is.null(sample_ids)) {
    sample_ids <- unique(as.character(cd[[sample_col]]))
  }
  set.seed(seed)
  out <- list()
  for (sample_id in sample_ids) {
    df <- cd[as.character(cd[[sample_col]]) == sample_id, , drop = FALSE]
    df <- df[!is.na(df[[x_col]]) & !is.na(df[[y_col]]), , drop = FALSE]
    if (nrow(df) == 0) next
    x_rng <- range(df[[x_col]], na.rm = TRUE)
    y_rng <- range(df[[y_col]], na.rm = TRUE)
    if (diff(x_rng) < roi_width_um || diff(y_rng) < roi_height_um) {
      warning("Sample ", sample_id, " is smaller than requested ROI size; skipping.")
      next
    }
    if (strategy == "random") {
      centers <- data.frame(
        x_center = stats::runif(grid_n * grid_n, x_rng[1] + roi_width_um / 2, x_rng[2] - roi_width_um / 2),
        y_center = stats::runif(grid_n * grid_n, y_rng[1] + roi_height_um / 2, y_rng[2] - roi_height_um / 2)
      )
    } else {
      centers <- expand.grid(
        x_center = seq(x_rng[1] + roi_width_um / 2, x_rng[2] - roi_width_um / 2, length.out = grid_n),
        y_center = seq(y_rng[1] + roi_height_um / 2, y_rng[2] - roi_height_um / 2, length.out = grid_n)
      )
    }
    candidates <- lapply(seq_len(nrow(centers)), function(i) {
      x_min <- centers$x_center[i] - roi_width_um / 2
      x_max <- centers$x_center[i] + roi_width_um / 2
      y_min <- centers$y_center[i] - roi_height_um / 2
      y_max <- centers$y_center[i] + roi_height_um / 2
      n_cells <- sum(df[[x_col]] >= x_min & df[[x_col]] <= x_max & df[[y_col]] >= y_min & df[[y_col]] <= y_max)
      data.frame(sample = sample_id, x_min = x_min, x_max = x_max, y_min = y_min, y_max = y_max,
                 x_center = centers$x_center[i], y_center = centers$y_center[i], n_cells = n_cells)
    })
    candidates <- do.call(rbind, candidates)
    candidates <- candidates[candidates$n_cells >= min_cells, , drop = FALSE]
    if (nrow(candidates) == 0) {
      warning("No ROI candidates met min_cells for sample ", sample_id, ".")
      next
    }
    if (strategy == "quantile_grid" && nrow(candidates) > n_rois) {
      candidates$density_rank <- rank(candidates$n_cells, ties.method = "random")
      targets <- stats::quantile(candidates$density_rank, probs = seq(0, 1, length.out = n_rois))
      keep <- unique(vapply(targets, function(z) which.min(abs(candidates$density_rank - z)), integer(1)))
      if (length(keep) < n_rois) {
        extra <- setdiff(order(candidates$n_cells, decreasing = TRUE), keep)
        keep <- c(keep, head(extra, n_rois - length(keep)))
      }
      candidates <- candidates[keep, , drop = FALSE]
    } else if (nrow(candidates) > n_rois) {
      candidates <- candidates[sample(seq_len(nrow(candidates)), n_rois), , drop = FALSE]
    }
    candidates$roi_id <- sprintf("%s_roi%02d", sample_id, seq_len(nrow(candidates)))
    out[[sample_id]] <- candidates[, c("roi_id", "sample", "x_min", "x_max", "y_min", "y_max", "x_center", "y_center", "n_cells")]
  }
  if (length(out) == 0) {
    return(data.frame(
      roi_id = character(), sample = character(), x_min = numeric(),
      x_max = numeric(), y_min = numeric(), y_max = numeric(),
      x_center = numeric(), y_center = numeric(), n_cells = integer()
    ))
  }
  result <- do.call(rbind, out)
  rownames(result) <- NULL
  result
}

#' Subset a Xenium CDS to cells falling inside selected ROIs
#'
#' @param cds A monocle3 cell_data_set.
#' @param rois ROI table from \code{select_xenium_rois()}.
#' @param sample_col,x_col,y_col Column names in \code{colData(cds)}.
#' @param roi_col Name of the added ROI membership column.
#' @return A subsetted cell_data_set with ROI membership in \code{colData}.
#' @export
subset_xenium_rois <- function(cds, rois, sample_col = "sample", x_col = "x_centroid", y_col = "y_centroid", roi_col = "roi_id") {
  cd <- as.data.frame(SummarizedExperiment::colData(cds))
  membership <- rep(NA_character_, nrow(cd))
  for (i in seq_len(nrow(rois))) {
    hit <- as.character(cd[[sample_col]]) == as.character(rois$sample[i]) &
      cd[[x_col]] >= rois$x_min[i] & cd[[x_col]] <= rois$x_max[i] &
      cd[[y_col]] >= rois$y_min[i] & cd[[y_col]] <= rois$y_max[i] &
      is.na(membership)
    membership[hit] <- as.character(rois$roi_id[i])
  }
  SummarizedExperiment::colData(cds)[[roi_col]] <- membership
  cds[, !is.na(membership)]
}

#' Plot Xenium ROI locations over cell centroids
#'
#' @param cds A monocle3 cell_data_set.
#' @param rois ROI table from \code{select_xenium_rois()}.
#' @param sample_id Optional sample to plot.
#' @param sample_col,x_col,y_col Column names in \code{colData(cds)}.
#' @param point_size,point_alpha Cell point appearance.
#' @return A ggplot object.
#' @export
plot_xenium_rois <- function(cds, rois, sample_id = NULL, sample_col = "sample", x_col = "x_centroid", y_col = "y_centroid", point_size = 0.05, point_alpha = 0.15) {
  cd <- as.data.frame(SummarizedExperiment::colData(cds))
  if (!is.null(sample_id)) {
    cd <- cd[as.character(cd[[sample_col]]) == sample_id, , drop = FALSE]
    rois <- rois[as.character(rois$sample) == sample_id, , drop = FALSE]
  }
  ggplot2::ggplot(cd, ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]])) +
    ggplot2::geom_point(size = point_size, alpha = point_alpha) +
    ggplot2::geom_rect(data = rois, ggplot2::aes(xmin = x_min, xmax = x_max, ymin = y_min, ymax = y_max, fill = roi_id),
                       inherit.aes = FALSE, alpha = 0.2, colour = "red", linewidth = 0.4) +
    ggplot2::geom_text(data = rois, ggplot2::aes(x = x_center, y = y_center, label = roi_id), inherit.aes = FALSE, size = 3) +
    ggplot2::coord_equal() +
    ggplot2::facet_wrap(stats::as.formula(paste("~", sample_col)), scales = "free") +
    ggplot2::theme_bw() +
    ggplot2::labs(x = paste0(x_col, " (um)"), y = paste0(y_col, " (um)"), fill = "ROI")
}

#' Summarize nearest neighboring cell types around focal cells
#'
#' @param cds A monocle3 cell_data_set, often already downsampled to ROIs.
#' @param focal_types Cell type(s) to query, e.g. \code{"POLC"}.
#' @param cell_type_col Broad annotation column in \code{colData(cds)}.
#' @param sample_col,x_col,y_col Column names in \code{colData(cds)}.
#' @param exclude_same_cell_type If TRUE, do not count neighbors with the same
#'   cell type as the focal cell.
#' @return A list with \code{nearest_cells} and \code{summary} data.frames.
#' @export
nearest_cell_type_summary <- function(cds, focal_types, cell_type_col, sample_col = "sample", x_col = "x_centroid", y_col = "y_centroid", exclude_same_cell_type = FALSE) {
  cd <- as.data.frame(SummarizedExperiment::colData(cds))
  required <- c(cell_type_col, sample_col, x_col, y_col)
  if (!all(required %in% colnames(cd))) stop("cds colData must contain: ", paste(required, collapse = ", "))
  cd$cell_id <- rownames(cd)
  cd <- cd[!is.na(cd[[cell_type_col]]) & !is.na(cd[[x_col]]) & !is.na(cd[[y_col]]), , drop = FALSE]
  hits <- list()
  idx <- 1L
  for (sample_id in unique(as.character(cd[[sample_col]]))) {
    df <- cd[as.character(cd[[sample_col]]) == sample_id, , drop = FALSE]
    focal <- df[df[[cell_type_col]] %in% focal_types, , drop = FALSE]
    for (i in seq_len(nrow(focal))) {
      pool <- df[df$cell_id != focal$cell_id[i], , drop = FALSE]
      if (exclude_same_cell_type) pool <- pool[pool[[cell_type_col]] != focal[[cell_type_col]][i], , drop = FALSE]
      if (nrow(pool) == 0) next
      d2 <- (pool[[x_col]] - focal[[x_col]][i])^2 + (pool[[y_col]] - focal[[y_col]][i])^2
      j <- which.min(d2)
      hits[[idx]] <- data.frame(sample = sample_id, focal_cell_id = focal$cell_id[i], focal_cell_type = focal[[cell_type_col]][i],
                                neighbor_cell_id = pool$cell_id[j], neighbor_cell_type = pool[[cell_type_col]][j], distance_um = sqrt(d2[j]))
      idx <- idx + 1L
    }
  }
  nearest <- if (length(hits)) do.call(rbind, hits) else data.frame()
  summary <- if (nrow(nearest)) {
    counts <- as.data.frame(table(nearest$focal_cell_type, nearest$neighbor_cell_type), stringsAsFactors = FALSE)
    colnames(counts) <- c("focal_cell_type", "neighbor_cell_type", "n")
    counts <- counts[counts$n > 0, , drop = FALSE]
    summed <- rowsum(counts$n, counts$focal_cell_type)
    totals <- stats::setNames(summed[, 1], rownames(summed))
    counts$proportion <- counts$n / totals[counts$focal_cell_type]
    counts[order(counts$focal_cell_type, -counts$n), , drop = FALSE]
  } else data.frame()
  list(nearest_cells = nearest, summary = summary)
}
