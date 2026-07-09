#' Read Xenium cell centroid coordinates and QC metrics
#'
#' Reads the cells.csv.gz file produced by Xenium Onboard Analysis,
#' which contains per-cell spatial centroid coordinates (in microns)
#' and per-cell QC metrics such as nucleus area and transcript counts.
#'
#' @param cells_csv Path to cells.csv.gz
#' @param barcode_suffix Optional cell ID suffix (e.g. "_1"), applied so
#'   cell IDs match the barcodes attached to a combined cell_data_set.
#' @param sample_id Optional sample identifier to attach as a "sample" column.
#'
#' @return A data.frame with spatial coordinates and QC metrics per cell
#' @export
read_xenium_cells <- function(
    cells_csv,
    barcode_suffix = NULL,
    sample_id      = NULL
) {

  if (!file.exists(cells_csv)) {
    stop("cells_csv does not exist: ", cells_csv)
  }

  dt <- data.table::fread(cells_csv)

  if (!"cell_id" %in% colnames(dt)) {
    stop("cells_csv must contain a 'cell_id' column")
  }

  if (!is.null(barcode_suffix)) {
    dt$cell_id <- paste0(dt$cell_id, barcode_suffix)
  }

  if (!is.null(sample_id)) {
    dt$sample <- sample_id
  }

  as.data.frame(dt)
}
