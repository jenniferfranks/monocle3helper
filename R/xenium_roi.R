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
