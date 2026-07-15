#' Analyze spatial proximity between Xenium cell types
#'
#' For each Xenium cell, finds the nearest neighboring cells within each sample
#' and summarizes which target cell types are closest to each source cell type.
#' Distances are calculated from centroid coordinates in microns, using the
#' \code{x_centroid} and \code{y_centroid} columns added by
#' \code{build_xenium_cds()} by default.
#'
#' @param cds A monocle3 cell_data_set with Xenium centroid coordinates and cell
#'   type annotations in \code{colData(cds)}.
#' @param cell_type_col Column in \code{colData(cds)} containing cell type labels.
#' @param sample_col Column in \code{colData(cds)} identifying independent Xenium
#'   samples. Nearest neighbors are only searched within the same sample. Set to
#'   \code{NULL} to analyze all cells together.
#' @param x_col,y_col Columns in \code{colData(cds)} containing centroid
#'   coordinates in microns.
#' @param k Number of nearest neighbors to keep per cell.
#' @param max_distance Optional maximum centroid distance in microns. Neighbor
#'   links farther than this value are discarded.
#' @param include_same_type Logical; if \code{FALSE}, neighbors with the same
#'   cell type as the source cell are excluded from the reported summaries.
#' @param chunk_size Number of source cells processed at once when calculating
#'   distances. Lower values use less memory.
#' @param return_neighbors Logical; if \code{TRUE}, include the per-cell nearest
#'   neighbor table in the returned list.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{summary}}{Per-sample source-to-target cell type proximity
#'   summary with neighbor counts and distance statistics.}
#'   \item{\code{overall_summary}}{The same statistics pooled across samples.}
#'   \item{\code{closest_cell_types}}{For each sample and source cell type, the
#'   target cell type with the smallest median nearest-neighbor distance.}
#'   \item{\code{neighbors}}{Returned only when \code{return_neighbors = TRUE};
#'   one row per retained source-neighbor link.}
#' }
#'
#' @export
analyze_xenium_cell_type_proximity <- function(
    cds,
    cell_type_col,
    sample_col = "sample",
    x_col = "x_centroid",
    y_col = "y_centroid",
    k = 1,
    max_distance = Inf,
    include_same_type = TRUE,
    chunk_size = 1000,
    return_neighbors = FALSE
) {
  if (!is.character(cell_type_col) || length(cell_type_col) != 1) {
    stop("cell_type_col must be a single column name.")
  }
  if (!is.null(sample_col) && (!is.character(sample_col) || length(sample_col) != 1)) {
    stop("sample_col must be NULL or a single column name.")
  }
  if (!is.numeric(k) || length(k) != 1 || is.na(k) || k < 1) {
    stop("k must be a positive integer.")
  }
  if (!is.numeric(chunk_size) || length(chunk_size) != 1 || is.na(chunk_size) || chunk_size < 1) {
    stop("chunk_size must be a positive integer.")
  }
  if (!is.numeric(max_distance) || length(max_distance) != 1 || is.na(max_distance) || max_distance < 0) {
    stop("max_distance must be a non-negative number or Inf.")
  }

  k <- as.integer(k)
  chunk_size <- as.integer(chunk_size)
  cd <- as.data.frame(SummarizedExperiment::colData(cds))
  required_cols <- c(cell_type_col, x_col, y_col)
  if (!is.null(sample_col)) {
    required_cols <- c(required_cols, sample_col)
  }
  missing_cols <- setdiff(required_cols, colnames(cd))
  if (length(missing_cols) > 0) {
    stop("Required columns missing from colData(cds): ", paste(missing_cols, collapse = ", "))
  }

  cd$cell_id <- rownames(cd)
  cd$.cell_type <- as.character(cd[[cell_type_col]])
  cd$.sample <- if (is.null(sample_col)) "all" else as.character(cd[[sample_col]])
  cd$.x <- as.numeric(cd[[x_col]])
  cd$.y <- as.numeric(cd[[y_col]])
  keep <- !is.na(cd$.cell_type) & nzchar(cd$.cell_type) &
    !is.na(cd$.sample) & nzchar(cd$.sample) &
    is.finite(cd$.x) & is.finite(cd$.y)
  dropped <- sum(!keep)
  if (dropped > 0) {
    warning(dropped, " cells with missing cell type, sample, or coordinates were ignored.")
  }
  cd <- cd[keep, , drop = FALSE]
  if (nrow(cd) < 2) {
    stop("At least two cells with complete metadata are required.")
  }

  neighbor_tables <- lapply(split(cd, cd$.sample), function(sample_cd) {
    find_nearest_xenium_neighbors(
      sample_cd = sample_cd,
      k = k,
      max_distance = max_distance,
      include_same_type = include_same_type,
      chunk_size = chunk_size
    )
  })
  neighbors <- dplyr::bind_rows(neighbor_tables)
  if (nrow(neighbors) == 0) {
    stop("No neighbor links were retained. Consider increasing max_distance or enabling include_same_type.")
  }

  summary <- summarize_xenium_proximity(neighbors, c("sample", "source_cell_type", "target_cell_type"))
  overall_summary <- summarize_xenium_proximity(neighbors, c("source_cell_type", "target_cell_type"))
  closest <- summary |>
    dplyr::group_by(sample, source_cell_type) |>
    dplyr::slice_min(median_distance, n = 1, with_ties = FALSE) |>
    dplyr::ungroup()

  result <- list(
    summary = summary,
    overall_summary = overall_summary,
    closest_cell_types = closest
  )
  if (return_neighbors) {
    result$neighbors <- neighbors
  }
  result
}

find_nearest_xenium_neighbors <- function(sample_cd, k, max_distance, include_same_type, chunk_size) {
  n <- nrow(sample_cd)
  if (n < 2) {
    return(data.frame())
  }
  xy <- as.matrix(sample_cd[, c(".x", ".y")])
  out <- list()
  out_i <- 1L
  for (start in seq.int(1L, n, by = chunk_size)) {
    end <- min(start + chunk_size - 1L, n)
    idx <- start:end
    dx <- outer(xy[idx, 1], xy[, 1], "-")
    dy <- outer(xy[idx, 2], xy[, 2], "-")
    dist <- sqrt(dx * dx + dy * dy)
    dist[cbind(seq_along(idx), idx)] <- Inf
    if (!include_same_type) {
      same_type <- outer(sample_cd$.cell_type[idx], sample_cd$.cell_type, "==")
      dist[same_type] <- Inf
    }
    for (row_i in seq_along(idx)) {
      finite <- which(is.finite(dist[row_i, ]) & dist[row_i, ] <= max_distance)
      if (length(finite) == 0) next
      ordered <- finite[order(dist[row_i, finite])]
      nearest <- head(ordered, k)
      source <- idx[row_i]
      out[[out_i]] <- data.frame(
        sample = sample_cd$.sample[source],
        source_cell_id = sample_cd$cell_id[source],
        source_cell_type = sample_cd$.cell_type[source],
        target_cell_id = sample_cd$cell_id[nearest],
        target_cell_type = sample_cd$.cell_type[nearest],
        neighbor_rank = seq_along(nearest),
        distance = dist[row_i, nearest],
        stringsAsFactors = FALSE
      )
      out_i <- out_i + 1L
    }
  }
  dplyr::bind_rows(out)
}

summarize_xenium_proximity <- function(neighbors, group_cols) {
  neighbors |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarise(
      n_source_cells = dplyr::n_distinct(source_cell_id),
      n_neighbor_links = dplyr::n(),
      mean_distance = mean(distance),
      median_distance = stats::median(distance),
      min_distance = min(distance),
      max_distance = max(distance),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::across(dplyr::all_of(group_cols)))
}
