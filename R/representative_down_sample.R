#' Downsample a cds to a target total number of cells, preserving rare cell types
#'
#' Downsamples a Monocle3 \code{cell_data_set} to (approximately) a target
#' total number of cells. Cell types are auto-detected from a metadata column
#' in \code{colData(cds)}. Cell types with fewer than or equal to
#' \code{min_cells} cells are kept in full (never downsampled). The remaining
#' \code{total_cells} budget (after reserving all cells from the small cell
#' types) is distributed across the larger cell types in proportion to their
#' size, and cells are randomly sampled (without replacement) within each
#' using \code{seed} for reproducibility.
#'
#' @param cds Monocle3 \code{cell_data_set}.
#' @param total_cells Target total number of cells in the downsampled cds.
#' @param min_cells Cell types with this many cells or fewer are kept in full
#'   and excluded from downsampling.
#' @param cell_type_col Column name in \code{colData(cds)} holding cell type
#'   labels. Defaults to \code{"cell_type"}.
#' @param seed Integer random seed used for reproducible sampling.
#'
#' @return A downsampled \code{cell_data_set} with (as close as possible to)
#'   \code{total_cells} cells.
#'
#' @importFrom SummarizedExperiment colData
#' @export
#'
#' @examples
#' \dontrun{
#' cds_small <- downsample_cds(
#'   cds,
#'   total_cells = 5000,
#'   min_cells = 100,
#'   cell_type_col = "cell_type",
#'   seed = 42
#' )
#' }
downsample_cds <- function(
    cds,
    total_cells,
    min_cells,
    cell_type_col = "cell_type",
    seed
) {
  if (!cell_type_col %in% colnames(colData(cds))) {
    stop(sprintf("'%s' not found in colData(cds)", cell_type_col))
  }
  n_total <- ncol(cds)
  if (total_cells >= n_total) {
    warning("total_cells >= number of cells in cds; returning cds unchanged")
    return(cds)
  }
  cell_types <- as.character(colData(cds)[[cell_type_col]])
  cell_ids <- colnames(cds)
  type_counts <- table(cell_types)
  small_types <- names(type_counts)[type_counts <= min_cells]
  large_types <- names(type_counts)[type_counts > min_cells]
  keep_ids <- cell_ids[cell_types %in% small_types]
  budget <- total_cells - length(keep_ids)
  if (budget < 0) {
    stop(
      "total_cells is smaller than the number of cells belonging to types ",
      "with <= min_cells cells; increase total_cells or lower min_cells"
    )
  }
  if (length(large_types) == 0 || budget == 0) {
    sampled_ids <- character(0)
  } else {
    large_counts <- type_counts[large_types]
    n_large <- sum(large_counts)
    # proportional allocation, capped at each type's own size
    alloc <- floor(budget * large_counts / n_large)
    alloc <- pmin(alloc, large_counts)
    # distribute any leftover cells (from flooring) one at a time to
    # the types with remaining room, largest-remaining-capacity first
    leftover <- budget - sum(alloc)
    if (leftover > 0) {
      room <- large_counts - alloc
      order_idx <- order(room, decreasing = TRUE)
      i <- 1
      while (leftover > 0 && any(room > 0)) {
        idx <- order_idx[i]
        if (room[idx] > 0) {
          alloc[idx] <- alloc[idx] + 1
          room[idx] <- room[idx] - 1
          leftover <- leftover - 1
        }
        i <- i %% length(order_idx) + 1
      }
    }
    set.seed(seed)
    sampled_ids <- unlist(lapply(large_types, function(ct) {
      ids_ct <- cell_ids[cell_types == ct]
      sample(ids_ct, size = alloc[[ct]])
    }), use.names = FALSE)
  }
  final_ids <- c(keep_ids, sampled_ids)
  cds[, final_ids]
}