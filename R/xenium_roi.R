#' Subset a Xenium CDS to rectangular regions of interest
#'
#' Creates a smaller \code{cell_data_set} containing only cells whose
#' Xenium centroid coordinates fall inside one or more rectangular regions of
#' interest (ROIs). Matching cells receive \code{roi_id} and
#' \code{roi_label} columns in \code{colData}. This is intended for creating
#' lightweight ROI-only CDS objects before running dimensionality reduction,
#' clustering, and exploratory analyses.
#'
#' @param cds A monocle3 cell_data_set with \code{x_centroid},
#'   \code{y_centroid}, and a sample column in \code{colData}.
#' @param roi_table A data.frame with one row per rectangular ROI. By default
#'   it must contain \code{sample}, \code{x_min}, \code{x_max},
#'   \code{y_min}, and \code{y_max}. Coordinates should be in the same
#'   micron coordinate system as Xenium \code{cells.csv.gz} centroids.
#' @param roi_id_col Optional column in \code{roi_table} containing ROI IDs.
#'   If \code{NULL}, IDs are generated as \code{roi_1}, \code{roi_2}, ...
#' @param roi_sample_col Column in \code{roi_table} containing sample IDs.
#' @param cds_sample_col Column in \code{colData(cds)} containing sample IDs.
#' @param x_min_col,x_max_col,y_min_col,y_max_col Columns in \code{roi_table}
#'   containing rectangle bounds in microns.
#' @param keep_unassigned Logical; if \code{FALSE}, return only cells inside
#'   an ROI. If \code{TRUE}, keep all cells and mark cells outside ROIs with
#'   \code{NA} ROI metadata.
#'
#' @return A monocle3 cell_data_set subset to ROI cells unless
#'   \code{keep_unassigned = TRUE}.
#'
#' @export
subset_xenium_rois <- function(
    cds,
    roi_table,
    roi_id_col     = NULL,
    roi_sample_col = "sample",
    cds_sample_col = "sample",
    x_min_col      = "x_min",
    x_max_col      = "x_max",
    y_min_col      = "y_min",
    y_max_col      = "y_max",
    keep_unassigned = FALSE
) {

  if (is.null(roi_table) || nrow(roi_table) == 0) {
    stop("roi_table must contain at least one ROI row.")
  }

  cd <- as.data.frame(SummarizedExperiment::colData(cds))
  required_cd <- c("x_centroid", "y_centroid", cds_sample_col)
  missing_cd <- setdiff(required_cd, colnames(cd))
  if (length(missing_cd) > 0) {
    stop(
      "cds is missing required colData columns: ",
      paste(missing_cd, collapse = ", ")
    )
  }

  required_roi <- c(roi_sample_col, x_min_col, x_max_col, y_min_col, y_max_col)
  missing_roi <- setdiff(required_roi, colnames(roi_table))
  if (length(missing_roi) > 0) {
    stop(
      "roi_table is missing required columns: ",
      paste(missing_roi, collapse = ", ")
    )
  }
  if (!is.null(roi_id_col) && !roi_id_col %in% colnames(roi_table)) {
    stop("roi_id_col is not present in roi_table: ", roi_id_col)
  }

  roi <- as.data.frame(roi_table)
  roi$.roi_id <- if (is.null(roi_id_col)) {
    paste0("roi_", seq_len(nrow(roi)))
  } else {
    as.character(roi[[roi_id_col]])
  }

  cd$.roi_id <- NA_character_

  for (i in seq_len(nrow(roi))) {
    sample_match <- as.character(cd[[cds_sample_col]]) ==
      as.character(roi[[roi_sample_col]][i])
    in_roi <- sample_match &
      !is.na(cd$x_centroid) & !is.na(cd$y_centroid) &
      cd$x_centroid >= roi[[x_min_col]][i] &
      cd$x_centroid <= roi[[x_max_col]][i] &
      cd$y_centroid >= roi[[y_min_col]][i] &
      cd$y_centroid <= roi[[y_max_col]][i]

    cd$.roi_id[in_roi & is.na(cd$.roi_id)] <- roi$.roi_id[i]
  }

  keep <- keep_unassigned | !is.na(cd$.roi_id)
  if (!any(keep)) {
    stop("No cells matched the ROI definitions in roi_table.")
  }

  cds_roi <- cds[, keep]
  matched_roi <- cd$.roi_id[keep]
  SummarizedExperiment::colData(cds_roi)$roi_id <- matched_roi
  SummarizedExperiment::colData(cds_roi)$roi_label <- matched_roi

  cds_roi
}

#' Derive rectangular Xenium ROI tables from existing data
#'
#' Automatically creates suggested rectangular ROIs per sample using the same
#' centroid-range strategy used by \code{build_xenium_images()} to define
#' cell-focused image extents. ROIs can be derived from an existing CDS,
#' directly from \code{sample_table$cells_csv}, or from the \code{images}
#' object returned by \code{build_xenium_images()}.
#'
#' When \code{n_rois_per_sample > 1}, the function suggests smaller ROIs by
#' placing rectangles inside each sample's observed cell/image extent. Set
#' \code{target_fraction} to the approximate fraction of each sample extent
#' that should be retained across all ROIs. By default ROIs do not overlap.
#'
#' @param cds Optional monocle3 cell_data_set with \code{x_centroid},
#'   \code{y_centroid}, and a sample column in \code{colData}.
#' @param sample_table Optional sample table produced by
#'   \code{discover_xenium_sample_table()}. If \code{cds} is \code{NULL},
#'   this is used to read each \code{cells_csv} and derive sample bounds
#'   without building a full CDS first.
#' @param images Optional object returned by \code{build_xenium_images()}.
#'   When supplied, \code{cells_xlim}/\code{cells_ylim} are converted from
#'   image pixels back into microns using \code{effective_pixel_size_um}.
#' @param sample_col Column in \code{colData(cds)} containing sample IDs.
#' @param n_rois_per_sample Number of ROIs to suggest for each sample.
#' @param target_fraction Approximate fraction of each sample's observed
#'   extent to retain across all ROIs. If \code{NULL}, defaults to 1 for a
#'   single ROI per sample and 0.1 for multiple ROIs per sample.
#' @param allow_overlap Logical; whether suggested ROIs may overlap. Defaults
#'   to \code{FALSE}. Non-overlapping ROIs may be smaller than requested if
#'   \code{target_fraction} and \code{n_rois_per_sample} cannot fit inside
#'   the sample extent without overlap.
#' @param max_cell_imbalance Maximum tolerated fractional imbalance in
#'   actual ROI cell counts within each sample before a warning is emitted.
#' @param roi_id_prefix Prefix for generated ROI IDs.
#'
#' @return A data.frame with columns \code{sample}, \code{roi_id},
#'   \code{x_min}, \code{x_max}, \code{y_min}, \code{y_max},
#'   \code{target_fraction}, \code{actual_fraction}, \code{n_cells}, and
#'   \code{n_cells_in_roi} when cell counts are available.
#'
#' @export
derive_xenium_roi_table <- function(
    cds = NULL,
    sample_table = NULL,
    images = NULL,
    sample_col = "sample",
    n_rois_per_sample = 1,
    target_fraction = NULL,
    allow_overlap = FALSE,
    max_cell_imbalance = 0.25,
    roi_id_prefix = "roi"
) {

  if (length(n_rois_per_sample) != 1 || n_rois_per_sample < 1) {
    stop("n_rois_per_sample must be a positive integer.")
  }
  n_rois_per_sample <- as.integer(n_rois_per_sample)

  if (is.null(target_fraction)) {
    target_fraction <- if (n_rois_per_sample == 1) 1 else 0.1
  }
  if (length(target_fraction) != 1 || target_fraction <= 0 || target_fraction > 1) {
    stop("target_fraction must be > 0 and <= 1.")
  }
  if (length(max_cell_imbalance) != 1 || max_cell_imbalance < 0) {
    stop("max_cell_imbalance must be a non-negative number.")
  }

  make_rois <- function(bounds, sample_id, sample_index, n_cells = NA_integer_) {
    full_width <- bounds$x_max - bounds$x_min
    full_height <- bounds$y_max - bounds$y_min
    if (!is.finite(full_width) || !is.finite(full_height) ||
        full_width <= 0 || full_height <= 0) {
      stop("Cannot derive ROI bounds for sample with invalid extent: ", sample_id)
    }

    full_area <- full_width * full_height
    roi_area <- full_area * target_fraction / n_rois_per_sample
    aspect <- full_width / full_height
    roi_width <- sqrt(roi_area * aspect)
    roi_height <- sqrt(roi_area / aspect)

    grid_cols <- ceiling(sqrt(n_rois_per_sample * aspect))
    grid_rows <- ceiling(n_rois_per_sample / grid_cols)
    cell_width <- full_width / grid_cols
    cell_height <- full_height / grid_rows

    if (!allow_overlap) {
      roi_width <- min(roi_width, cell_width)
      roi_height <- min(roi_height, cell_height)
    } else {
      roi_width <- min(roi_width, full_width)
      roi_height <- min(roi_height, full_height)
    }

    rows <- vector("list", n_rois_per_sample)
    for (j in seq_len(n_rois_per_sample)) {
      grid_col <- ((j - 1) %% grid_cols) + 1
      grid_row <- floor((j - 1) / grid_cols) + 1
      center_x <- bounds$x_min + (grid_col - 0.5) * cell_width
      center_y <- bounds$y_min + (grid_row - 0.5) * cell_height

      x_min <- max(bounds$x_min, center_x - roi_width / 2)
      x_max <- min(bounds$x_max, center_x + roi_width / 2)
      y_min <- max(bounds$y_min, center_y - roi_height / 2)
      y_max <- min(bounds$y_max, center_y + roi_height / 2)

      rows[[j]] <- data.frame(
        sample = sample_id,
        roi_id = paste0(roi_id_prefix, "_", sample_index, "_", j),
        x_min = x_min,
        x_max = x_max,
        y_min = y_min,
        y_max = y_max,
        target_fraction = target_fraction,
        actual_fraction = ((x_max - x_min) * (y_max - y_min)) / full_area,
        n_cells = n_cells,
        stringsAsFactors = FALSE
      )
    }
    do.call(rbind, rows)
  }

  make_cell_balanced_rois <- function(df_samp, sample_id, sample_index) {
    df_samp <- df_samp[
      is.finite(df_samp$x_centroid) & is.finite(df_samp$y_centroid),
      ,
      drop = FALSE
    ]
    n_sample_cells <- nrow(df_samp)
    if (n_sample_cells < n_rois_per_sample) {
      stop(
        "Sample ", sample_id, " has ", n_sample_cells,
        " cells with finite centroids, fewer than n_rois_per_sample = ",
        n_rois_per_sample, ". Reduce n_rois_per_sample."
      )
    }

    total_target_cells <- max(
      n_rois_per_sample,
      min(n_sample_cells, round(n_sample_cells * target_fraction))
    )
    base_cells <- total_target_cells %/% n_rois_per_sample
    extra_cells <- total_target_cells %% n_rois_per_sample
    target_counts <- rep(base_cells, n_rois_per_sample)
    if (extra_cells > 0) {
      target_counts[seq_len(extra_cells)] <- target_counts[seq_len(extra_cells)] + 1
    }

    ord <- order(df_samp$x_centroid, df_samp$y_centroid)
    df_sorted <- df_samp[ord, , drop = FALSE]
    bin_breaks <- round(seq(0, n_sample_cells, length.out = n_rois_per_sample + 1))

    bounds <- list(
      x_min = min(df_samp$x_centroid, na.rm = TRUE),
      x_max = max(df_samp$x_centroid, na.rm = TRUE),
      y_min = min(df_samp$y_centroid, na.rm = TRUE),
      y_max = max(df_samp$y_centroid, na.rm = TRUE)
    )
    full_area <- (bounds$x_max - bounds$x_min) * (bounds$y_max - bounds$y_min)

    rows <- vector("list", n_rois_per_sample)
    actual_counts <- integer(n_rois_per_sample)

    for (j in seq_len(n_rois_per_sample)) {
      bin_start <- bin_breaks[j] + 1
      bin_end <- bin_breaks[j + 1]
      bin_index <- bin_start:bin_end
      bin_index <- bin_index[bin_index >= 1 & bin_index <= n_sample_cells]
      if (length(bin_index) == 0) {
        stop("Internal ROI binning produced an empty bin for sample: ", sample_id)
      }

      if (allow_overlap) {
        target_n <- target_counts[j]
        center_pos <- round(stats::quantile(
          seq_len(n_sample_cells),
          probs = j / (n_rois_per_sample + 1),
          names = FALSE
        ))
        start_pos <- max(1, center_pos - floor((target_n - 1) / 2))
        end_pos <- start_pos + target_n - 1
        if (end_pos > n_sample_cells) {
          end_pos <- n_sample_cells
          start_pos <- end_pos - target_n + 1
        }
      } else {
        target_n <- min(target_counts[j], length(bin_index))
        center_pos <- floor(stats::median(bin_index))
        start_pos <- max(min(bin_index), center_pos - floor((target_n - 1) / 2))
        end_pos <- start_pos + target_n - 1
        if (end_pos > max(bin_index)) {
          end_pos <- max(bin_index)
          start_pos <- end_pos - target_n + 1
        }
      }
      selected <- df_sorted[start_pos:end_pos, , drop = FALSE]

      x_min <- min(selected$x_centroid, na.rm = TRUE)
      x_max <- max(selected$x_centroid, na.rm = TRUE)
      y_min <- min(selected$y_centroid, na.rm = TRUE)
      y_max <- max(selected$y_centroid, na.rm = TRUE)

      # Give one-cell-wide ROIs a tiny data-derived width/height so downstream
      # inclusive bound checks retain the selected cell without extending
      # outside the sample extent.
      if (x_min == x_max) {
        eps_x <- max((bounds$x_max - bounds$x_min) * 1e-6, .Machine$double.eps)
        x_min <- max(bounds$x_min, x_min - eps_x)
        x_max <- min(bounds$x_max, x_max + eps_x)
      }
      if (y_min == y_max) {
        eps_y <- max((bounds$y_max - bounds$y_min) * 1e-6, .Machine$double.eps)
        y_min <- max(bounds$y_min, y_min - eps_y)
        y_max <- min(bounds$y_max, y_max + eps_y)
      }

      actual_counts[j] <- sum(
        df_samp$x_centroid >= x_min & df_samp$x_centroid <= x_max &
          df_samp$y_centroid >= y_min & df_samp$y_centroid <= y_max,
        na.rm = TRUE
      )

      rows[[j]] <- data.frame(
        sample = sample_id,
        roi_id = paste0(roi_id_prefix, "_", sample_index, "_", j),
        x_min = x_min,
        x_max = x_max,
        y_min = y_min,
        y_max = y_max,
        target_fraction = target_fraction,
        actual_fraction = ((x_max - x_min) * (y_max - y_min)) / full_area,
        n_cells = n_sample_cells,
        n_cells_in_roi = actual_counts[j],
        stringsAsFactors = FALSE
      )
    }

    if (any(actual_counts == 0)) {
      stop("At least one suggested ROI contains zero cells for sample: ", sample_id)
    }

    mean_count <- mean(actual_counts)
    if (mean_count > 0 && any(abs(actual_counts - mean_count) / mean_count > max_cell_imbalance)) {
      warning(
        "Suggested ROI cell counts vary by more than ",
        max_cell_imbalance * 100,
        "% for sample ", sample_id,
        ". Inspect n_cells_in_roi and adjust n_rois_per_sample or target_fraction."
      )
    }

    do.call(rbind, rows)
  }

  if (!is.null(images)) {
    if (is.null(images$images) || is.null(images$effective_pixel_size_um)) {
      stop("images must be an object returned by build_xenium_images().")
    }

    effective_pixel_size <- images$effective_pixel_size_um
    image_names <- names(images$images)
    roi_rows <- lapply(seq_along(image_names), function(i) {
      sample_id <- image_names[i]
      img <- images$images[[sample_id]]
      if (is.null(img$cells_xlim) || is.null(img$cells_ylim)) {
        stop("images$images entries must contain cells_xlim and cells_ylim.")
      }
      bounds <- list(
        x_min = min(img$cells_xlim) * effective_pixel_size,
        x_max = max(img$cells_xlim) * effective_pixel_size,
        y_min = min(img$cells_ylim) * effective_pixel_size,
        y_max = max(img$cells_ylim) * effective_pixel_size
      )
      make_rois(bounds, sample_id, i)
    })

    return(do.call(rbind, roi_rows))
  }

  if (!is.null(cds)) {
    cd <- as.data.frame(SummarizedExperiment::colData(cds))
    required_cd <- c("x_centroid", "y_centroid", sample_col)
    missing_cd <- setdiff(required_cd, colnames(cd))
    if (length(missing_cd) > 0) {
      stop(
        "cds is missing required colData columns: ",
        paste(missing_cd, collapse = ", ")
      )
    }
    cd$.sample <- as.character(cd[[sample_col]])
  } else if (!is.null(sample_table)) {
    required_sample <- c("sample_id", "cells_csv", "barcode_suffix")
    missing_sample <- setdiff(required_sample, colnames(sample_table))
    if (length(missing_sample) > 0) {
      stop(
        "sample_table is missing required columns: ",
        paste(missing_sample, collapse = ", ")
      )
    }

    cell_meta <- vector("list", nrow(sample_table))
    for (i in seq_len(nrow(sample_table))) {
      cell_meta[[i]] <- read_xenium_cells(
        cells_csv = sample_table$cells_csv[i],
        barcode_suffix = sample_table$barcode_suffix[i],
        sample_id = sample_table$sample_id[i]
      )
    }
    cd <- dplyr::bind_rows(cell_meta)
    cd$.sample <- as.character(cd$sample)
  } else {
    stop("Provide cds, sample_table, or images to derive ROI bounds.")
  }

  samples <- unique(cd$.sample)
  roi_rows <- lapply(seq_along(samples), function(i) {
    sample_id <- samples[i]
    df_samp <- cd[cd$.sample == sample_id, , drop = FALSE]
    make_cell_balanced_rois(df_samp, sample_id, i)
  })

  roi_table <- do.call(rbind, roi_rows)

  if (!is.null(sample_table)) {
    sample_cols <- intersect(
      c("sample_id", "sample_label", "sample_group"),
      colnames(sample_table)
    )
    sample_info <- unique(as.data.frame(sample_table[sample_cols]))
    roi_table <- dplyr::left_join(
      roi_table,
      sample_info,
      by = c("sample" = "sample_id")
    )
  }

  roi_table
}
