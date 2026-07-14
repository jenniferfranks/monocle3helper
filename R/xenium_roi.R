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

#' Derive a rectangular Xenium ROI table from existing data
#'
#' Automatically creates one rectangular ROI per sample using the same
#' centroid-range strategy used by \code{build_xenium_images()} to define
#' cell-focused image extents. ROIs can be derived from an existing CDS,
#' directly from \code{sample_table$cells_csv}, or from the \code{images}
#' object returned by \code{build_xenium_images()}.
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
#' @param buffer_um Extra micron buffer to add around centroid-derived ROIs.
#'   Ignored when \code{images} is supplied because image extents already
#'   include the image buffer used by \code{build_xenium_images()}.
#' @param roi_id_prefix Prefix for generated ROI IDs.
#'
#' @return A data.frame with columns \code{sample}, \code{roi_id},
#'   \code{x_min}, \code{x_max}, \code{y_min}, \code{y_max}, and
#'   \code{n_cells} when cell counts are available.
#'
#' @export
derive_xenium_roi_table <- function(
    cds = NULL,
    sample_table = NULL,
    images = NULL,
    sample_col = "sample",
    buffer_um = 0,
    roi_id_prefix = "roi"
) {

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
      data.frame(
        sample = sample_id,
        roi_id = paste0(roi_id_prefix, "_", i),
        x_min = min(img$cells_xlim) * effective_pixel_size,
        x_max = max(img$cells_xlim) * effective_pixel_size,
        y_min = min(img$cells_ylim) * effective_pixel_size,
        y_max = max(img$cells_ylim) * effective_pixel_size,
        n_cells = NA_integer_,
        stringsAsFactors = FALSE
      )
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
    data.frame(
      sample = sample_id,
      roi_id = paste0(roi_id_prefix, "_", i),
      x_min = min(df_samp$x_centroid, na.rm = TRUE) - buffer_um,
      x_max = max(df_samp$x_centroid, na.rm = TRUE) + buffer_um,
      y_min = min(df_samp$y_centroid, na.rm = TRUE) - buffer_um,
      y_max = max(df_samp$y_centroid, na.rm = TRUE) + buffer_um,
      n_cells = nrow(df_samp),
      stringsAsFactors = FALSE
    )
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
