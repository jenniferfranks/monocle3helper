#' Discover Xenium sample metadata on disk
#'
#' Automatically scans a root directory for 10x Genomics Xenium
#' Onboard Analysis outputs and constructs a per-sample metadata
#' table describing file locations and sample annotations.
#'
#' Each sample directory is expected to contain:
#' \itemize{
#'   \item cell_feature_matrix.h5
#'   \item cells.csv.gz
#'   \item morphology_focus/ (directory of per-channel OME-TIFF images)
#' }
#'
#' @param root_dir Directory containing per-sample Xenium output folders.
#' @param sample_dirs Optional character vector of sample subdirectories
#'   (relative to \code{root_dir} or full paths). If NULL, all first-level
#'   subdirectories of \code{root_dir} are used.
#' @param sample_ids Optional character vector of full sample IDs.
#'   Defaults to directory names.
#' @param sample_labels Optional vector of colloquial sample labels.
#'   Defaults to \code{sample_ids}.
#' @param sample_groups Optional vector of grouping labels.
#'   Defaults to \code{sample_labels}.
#' @param barcode_suffix Optional vector of cell ID suffixes (for example "_1").
#'   Defaults to sequential suffixes.
#'
#' @return A tibble with one row per sample and columns:
#' \describe{
#'   \item{sample_id}{Full sample identifier}
#'   \item{sample_label}{Colloquial sample label}
#'   \item{sample_group}{Grouping label}
#'   \item{h5_file}{Path to cell_feature_matrix.h5}
#'   \item{cells_csv}{Path to cells.csv.gz}
#'   \item{morphology_dir}{Path to the morphology_focus directory}
#'   \item{barcode_suffix}{Cell ID suffix used for this sample}
#' }
#'
#' @export
discover_xenium_sample_table <- function(
    root_dir,
    sample_dirs    = NULL,
    sample_ids     = NULL,
    sample_labels  = NULL,
    sample_groups  = NULL,
    barcode_suffix = NULL
) {

  if (!dir.exists(root_dir)) {
    stop("root_dir does not exist: ", root_dir)
  }

  ## ------------------------------------------------------------
  ## Determine sample directories
  ## ------------------------------------------------------------
  if (is.null(sample_dirs)) {
    sample_paths <- list.dirs(
      path       = root_dir,
      full.names = TRUE,
      recursive  = FALSE
    )
  } else {
    sample_paths <- ifelse(
      dir.exists(sample_dirs),
      sample_dirs,
      file.path(root_dir, sample_dirs)
    )
  }

  if (length(sample_paths) == 0) {
    stop("No sample directories found under: ", root_dir)
  }

  n <- length(sample_paths)

  ## ------------------------------------------------------------
  ## Sample identifiers
  ## ------------------------------------------------------------
  if (is.null(sample_ids)) {
    sample_ids <- basename(sample_paths)
  } else if (length(sample_ids) != n) {
    stop("sample_ids must have length equal to number of samples")
  }

  if (is.null(sample_labels)) {
    sample_labels <- sample_ids
  } else if (length(sample_labels) != n) {
    stop("sample_labels must have length equal to number of samples")
  }

  if (is.null(sample_groups)) {
    sample_groups <- sample_labels
  } else if (length(sample_groups) != n) {
    stop("sample_groups must have length equal to number of samples")
  }

  if (is.null(barcode_suffix)) {
    barcode_suffix <- paste0("_", seq_len(n))
  } else if (length(barcode_suffix) != n) {
    stop("barcode_suffix must have length equal to number of samples")
  }

  ## ------------------------------------------------------------
  ## Discover per-sample files
  ## ------------------------------------------------------------
  h5_files       <- character(n)
  cells_files    <- character(n)
  morphology_dirs <- character(n)

  for (i in seq_len(n)) {
    sp <- sample_paths[i]

    h5_files[i] <- list.files(
      sp,
      pattern    = "cell_feature_matrix\\.h5$",
      full.names = TRUE,
      recursive  = TRUE
    )[1]

    cells_files[i] <- list.files(
      sp,
      pattern    = "cells\\.csv\\.gz$",
      full.names = TRUE,
      recursive  = TRUE
    )[1]

    morph_dir <- list.dirs(
      sp,
      full.names = TRUE,
      recursive  = TRUE
    )
    morph_dir <- morph_dir[basename(morph_dir) == "morphology_focus"]
    morphology_dirs[i] <- if (length(morph_dir) > 0) morph_dir[1] else NA_character_
  }

  ## ------------------------------------------------------------
  ## Assemble result table
  ## ------------------------------------------------------------
  tibble::tibble(
    sample_id      = sample_ids,
    sample_label   = sample_labels,
    sample_group   = sample_groups,
    h5_file        = h5_files,
    cells_csv      = cells_files,
    morphology_dir = morphology_dirs,
    barcode_suffix = barcode_suffix
  )
}
