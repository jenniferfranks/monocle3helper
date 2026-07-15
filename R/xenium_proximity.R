#' Analyze spatial proximity between Xenium cell types
#'
#' Finds nearest neighbors independently within each stratum (for example,
#' sample or sample + ROI), reports directed focal cell type to neighbor cell
#' type relationships, summarizes spatial distances, and adjusts associations
#' for neighbor cell type abundance using enrichment ratios.
#'
#' @param cds A monocle3 cell_data_set with Xenium centroid coordinates and cell
#'   type annotations in \code{colData(cds)}.
#' @param cell_type_col Column in \code{colData(cds)} containing cell type labels.
#' @param strata_cols Character vector of \code{colData(cds)} columns that define
#'   independent neighborhoods. Nearest neighbors are never searched across
#'   strata. Defaults to \code{"sample"}.
#' @param k Number of nearest neighbors to retain for each focal cell.
#' @param max_distance Maximum centroid-to-centroid distance in microns for
#'   retained neighbors. Use \code{Inf} to keep the nearest \code{k} neighbors
#'   regardless of distance.
#' @param x_col,y_col Columns in \code{colData(cds)} containing centroid
#'   coordinates in microns.
#' @param roi_col Optional ROI column to append to \code{strata_cols}. This is a
#'   convenience for analyzing neighborhoods independently within ROIs created
#'   by \code{subset_xenium_rois()}.
#' @param exclude_same_type Logical; if \code{TRUE}, same-cell-type neighbor
#'   interactions are excluded before summaries and enrichment are calculated.
#' @param return_cell_neighbors Logical; if \code{TRUE}, include one row per
#'   retained focal cell to neighbor cell link.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{pairwise}}{Directed focal type to neighbor type summaries by
#'   stratum, including observed links, expected links, enrichment, mean
#'   distance, and median distance.}
#'   \item{\code{by_stratum}}{Cell type abundance and retained-neighbor counts
#'   for each stratum.}
#'   \item{\code{cell_neighbors}}{Per-cell neighbor links, returned when
#'   \code{return_cell_neighbors = TRUE}.}
#' }
#'
#' @export
analyze_xenium_cell_proximity <- function(
    cds,
    cell_type_col = "cell_type",
    strata_cols = "sample",
    k = 10,
    max_distance = 50,
    x_col = "x_centroid",
    y_col = "y_centroid",
    roi_col = NULL,
    exclude_same_type = FALSE,
    return_cell_neighbors = TRUE
) {
  if (!requireNamespace("RANN", quietly = TRUE)) {
    stop("Package 'RANN' is required for analyze_xenium_cell_proximity().")
  }
  validate_xenium_proximity_args(
    cell_type_col = cell_type_col,
    strata_cols = strata_cols,
    k = k,
    max_distance = max_distance,
    x_col = x_col,
    y_col = y_col,
    roi_col = roi_col
  )

  if (!is.null(roi_col)) {
    strata_cols <- unique(c(strata_cols, roi_col))
  }

  cd <- as.data.frame(SummarizedExperiment::colData(cds))
  required_cols <- unique(c(cell_type_col, strata_cols, x_col, y_col))
  missing_cols <- setdiff(required_cols, colnames(cd))
  if (length(missing_cols) > 0) {
    stop("Required columns missing from colData(cds): ", paste(missing_cols, collapse = ", "))
  }

  cd$cell_id <- rownames(cd)
  cd$.cell_type <- as.character(cd[[cell_type_col]])
  cd$.x <- as.numeric(cd[[x_col]])
  cd$.y <- as.numeric(cd[[y_col]])
  cd$.stratum <- make_xenium_proximity_stratum(cd, strata_cols)

  keep <- !is.na(cd$.cell_type) & nzchar(cd$.cell_type) &
    !is.na(cd$.stratum) & nzchar(cd$.stratum) &
    is.finite(cd$.x) & is.finite(cd$.y)
  dropped <- sum(!keep)
  if (dropped > 0) {
    warning(dropped, " cells with missing cell type, stratum, or coordinates were ignored.")
  }
  cd <- cd[keep, , drop = FALSE]
  if (nrow(cd) < 2) {
    stop("At least two cells with complete metadata are required.")
  }

  split_cd <- split(cd, cd$.stratum)
  neighbor_tables <- lapply(split_cd, find_xenium_cell_neighbors,
    strata_cols = strata_cols,
    k = as.integer(k),
    max_distance = max_distance,
    exclude_same_type = exclude_same_type
  )
  cell_neighbors <- dplyr::bind_rows(neighbor_tables)
  if (nrow(cell_neighbors) == 0) {
    stop("No neighbor links were retained. Increase max_distance or disable exclude_same_type.")
  }

  by_stratum <- summarize_xenium_proximity_strata(cd, cell_neighbors, strata_cols)
  pairwise <- summarize_xenium_pairwise_proximity(
    cell_neighbors = cell_neighbors,
    by_stratum = by_stratum,
    strata_cols = strata_cols,
    exclude_same_type = exclude_same_type
  )

  result <- list(pairwise = pairwise, by_stratum = by_stratum)
  if (return_cell_neighbors) {
    result$cell_neighbors <- cell_neighbors
  }
  attr(result, "strata_cols") <- strata_cols
  class(result) <- c("xenium_cell_proximity", class(result))
  result
}

#' @export
analyze_xenium_cell_type_proximity <- function(
    cds,
    cell_type_col,
    sample_col = "sample",
    x_col = "x_centroid",
    y_col = "y_centroid",
    k = 10,
    max_distance = 50,
    include_same_type = TRUE,
    chunk_size = 1000,
    return_neighbors = TRUE
) {
  warning(
    "analyze_xenium_cell_type_proximity() is deprecated; use ",
    "analyze_xenium_cell_proximity() instead.",
    call. = FALSE
  )
  analyze_xenium_cell_proximity(
    cds = cds,
    cell_type_col = cell_type_col,
    strata_cols = sample_col,
    k = k,
    max_distance = max_distance,
    x_col = x_col,
    y_col = y_col,
    exclude_same_type = !include_same_type,
    return_cell_neighbors = return_neighbors
  )
}

#' Plot Xenium cell type proximity summaries
#'
#' @param proximity Object returned by \code{analyze_xenium_cell_proximity()}.
#' @param value Pairwise metric to plot. Common choices are \code{"enrichment"}
#'   and \code{"median_distance"}.
#' @param facet Logical; if \code{TRUE}, facet the heatmap by stratum.
#'
#' @return A ggplot heatmap.
#'
#' @export
plot_xenium_cell_proximity <- function(
    proximity,
    value = c("enrichment", "median_distance"),
    facet = TRUE
) {
  if (!inherits(proximity, "xenium_cell_proximity") || is.null(proximity$pairwise)) {
    stop("proximity must be an object returned by analyze_xenium_cell_proximity().")
  }
  value <- match.arg(value)
  pairwise <- proximity$pairwise
  if (!value %in% colnames(pairwise)) {
    stop("value is not present in proximity$pairwise: ", value)
  }

  p <- ggplot2::ggplot(
    pairwise,
    ggplot2::aes_string(
      x = "neighbor_cell_type",
      y = "focal_cell_type",
      fill = value
    )
  ) +
    ggplot2::geom_tile() +
    ggplot2::labs(
      x = "Neighbor cell type",
      y = "Focal cell type",
      fill = value
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  if (facet && "stratum" %in% colnames(pairwise)) {
    p <- p + ggplot2::facet_wrap(stats::as.formula("~ stratum"))
  }
  p
}

#' Compare Xenium cell proximity statistics between strata
#'
#' Computes focal cell type to neighbor cell type differences between sample or
#' ROI strata from an object returned by \code{analyze_xenium_cell_proximity()}.
#' Use this after preserving sample and/or ROI columns in \code{strata_cols}.
#'
#' @param proximity Object returned by \code{analyze_xenium_cell_proximity()}.
#' @param stratum_col Column in \code{proximity$pairwise} to compare, such as
#'   \code{"sample"} or \code{"roi_id"}.
#' @param reference Reference stratum level. If \code{NULL}, the first sorted
#'   level is used.
#' @param value Pairwise statistic to compare.
#' @param within_cols Optional columns that define matched contexts for the
#'   comparison. For example, use \code{"roi_id"} when comparing samples within
#'   each ROI.
#'
#' @return A data.frame containing reference values, comparison values,
#'   differences, and log2 ratios for each directed cell-type pair.
#'
#' @export
compare_xenium_cell_proximity <- function(
    proximity,
    stratum_col = "sample",
    reference = NULL,
    value = c("enrichment", "median_distance", "observed_fraction"),
    within_cols = NULL
) {
  if (!inherits(proximity, "xenium_cell_proximity") || is.null(proximity$pairwise)) {
    stop("proximity must be an object returned by analyze_xenium_cell_proximity().")
  }
  value <- match.arg(value)
  pairwise <- proximity$pairwise
  required_cols <- c(stratum_col, within_cols, "focal_cell_type", "neighbor_cell_type", value)
  missing_cols <- setdiff(required_cols, colnames(pairwise))
  if (length(missing_cols) > 0) {
    stop("Columns missing from proximity$pairwise: ", paste(missing_cols, collapse = ", "))
  }

  levels <- sort(unique(as.character(pairwise[[stratum_col]])))
  if (length(levels) < 2) {
    stop("At least two levels in stratum_col are required for comparison.")
  }
  if (is.null(reference)) {
    reference <- levels[1]
  }
  if (!reference %in% levels) {
    stop("reference is not present in ", stratum_col, ": ", reference)
  }

  key_cols <- c(within_cols, "focal_cell_type", "neighbor_cell_type")
  ref <- pairwise[pairwise[[stratum_col]] == reference, c(key_cols, value), drop = FALSE]
  names(ref)[names(ref) == value] <- "reference_value"

  comp <- pairwise[pairwise[[stratum_col]] != reference, c(stratum_col, key_cols, value), drop = FALSE]
  names(comp)[names(comp) == value] <- "comparison_value"

  out <- dplyr::left_join(comp, ref, by = key_cols)
  out$reference <- reference
  out$value <- value
  out$difference <- out$comparison_value - out$reference_value
  out$log2_ratio <- log2((out$comparison_value + .Machine$double.eps) /
    (out$reference_value + .Machine$double.eps))
  out[, c(stratum_col, "reference", within_cols, "focal_cell_type", "neighbor_cell_type",
    "value", "comparison_value", "reference_value", "difference", "log2_ratio"), drop = FALSE]
}

validate_xenium_proximity_args <- function(cell_type_col, strata_cols, k, max_distance, x_col, y_col, roi_col) {
  if (!is.character(cell_type_col) || length(cell_type_col) != 1) {
    stop("cell_type_col must be a single column name.")
  }
  if (!is.character(strata_cols) || length(strata_cols) < 1 || anyNA(strata_cols)) {
    stop("strata_cols must be a non-empty character vector.")
  }
  if (!is.numeric(k) || length(k) != 1 || is.na(k) || k < 1) {
    stop("k must be a positive integer.")
  }
  if (!is.numeric(max_distance) || length(max_distance) != 1 || is.na(max_distance) || max_distance < 0) {
    stop("max_distance must be a non-negative number or Inf.")
  }
  if (!is.character(x_col) || length(x_col) != 1 || !is.character(y_col) || length(y_col) != 1) {
    stop("x_col and y_col must each be a single column name.")
  }
  if (!is.null(roi_col) && (!is.character(roi_col) || length(roi_col) != 1)) {
    stop("roi_col must be NULL or a single column name.")
  }
}

make_xenium_proximity_stratum <- function(cd, strata_cols) {
  strata_df <- cd[, strata_cols, drop = FALSE]
  complete <- stats::complete.cases(strata_df)
  stratum <- rep(NA_character_, nrow(cd))
  stratum[complete] <- apply(strata_df[complete, , drop = FALSE], 1, paste, collapse = " | ")
  stratum
}

find_xenium_cell_neighbors <- function(sample_cd, strata_cols, k, max_distance, exclude_same_type) {
  n <- nrow(sample_cd)
  if (n < 2) {
    return(data.frame())
  }

  xy <- as.matrix(sample_cd[, c(".x", ".y")])
  k_search <- if (exclude_same_type) n else min(n, k + 1L)
  nn <- RANN::nn2(data = xy, query = xy, k = k_search, searchtype = "standard")

  out <- vector("list", n)
  for (i in seq_len(n)) {
    neighbor_idx <- nn$nn.idx[i, ]
    distances <- nn$nn.dists[i, ]
    keep <- neighbor_idx != i & is.finite(distances) & distances <= max_distance
    if (exclude_same_type) {
      keep <- keep & sample_cd$.cell_type[neighbor_idx] != sample_cd$.cell_type[i]
    }
    neighbor_idx <- neighbor_idx[keep]
    distances <- distances[keep]
    if (length(neighbor_idx) == 0) next

    ord <- order(distances)
    neighbor_idx <- head(neighbor_idx[ord], k)
    distances <- head(distances[ord], k)
    stratum_values <- sample_cd[rep(i, length(neighbor_idx)), strata_cols, drop = FALSE]
    out[[i]] <- data.frame(
      stratum_values,
      stratum = sample_cd$.stratum[i],
      focal_cell_id = sample_cd$cell_id[i],
      focal_cell_type = sample_cd$.cell_type[i],
      neighbor_cell_id = sample_cd$cell_id[neighbor_idx],
      neighbor_cell_type = sample_cd$.cell_type[neighbor_idx],
      neighbor_rank = seq_along(neighbor_idx),
      distance = distances,
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }
  dplyr::bind_rows(out)
}

summarize_xenium_proximity_strata <- function(cd, cell_neighbors, strata_cols) {
  abundance <- cd |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(strata_cols, "stratum"))), .cell_type) |>
    dplyr::summarise(n_cells = dplyr::n(), .groups = "drop_last") |>
    dplyr::mutate(stratum_total_cells = sum(n_cells), abundance_fraction = n_cells / stratum_total_cells) |>
    dplyr::ungroup() |>
    dplyr::rename(cell_type = .cell_type)

  retained <- cell_neighbors |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(strata_cols, "stratum"))), focal_cell_type) |>
    dplyr::summarise(retained_neighbor_links = dplyr::n(), .groups = "drop") |>
    dplyr::rename(cell_type = focal_cell_type)

  dplyr::left_join(abundance, retained, by = c(strata_cols, "stratum", "cell_type")) |>
    dplyr::mutate(retained_neighbor_links = tidyr_coalesce_zero(retained_neighbor_links))
}

summarize_xenium_pairwise_proximity <- function(cell_neighbors, by_stratum, strata_cols, exclude_same_type) {
  pairwise <- cell_neighbors |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(strata_cols, "stratum", "focal_cell_type", "neighbor_cell_type")))) |>
    dplyr::summarise(
      observed_links = dplyr::n(),
      n_focal_cells_with_neighbor = dplyr::n_distinct(focal_cell_id),
      mean_distance = mean(distance),
      median_distance = stats::median(distance),
      .groups = "drop"
    )

  focal_totals <- cell_neighbors |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(strata_cols, "stratum", "focal_cell_type")))) |>
    dplyr::summarise(total_links_from_focal_type = dplyr::n(), .groups = "drop")

  abundance <- by_stratum |>
    dplyr::select(dplyr::all_of(c(strata_cols, "stratum", "cell_type", "n_cells", "stratum_total_cells"))) |>
    dplyr::rename(neighbor_cell_type = cell_type, neighbor_type_cells = n_cells)

  focal_abundance <- by_stratum |>
    dplyr::select(dplyr::all_of(c(strata_cols, "stratum", "cell_type", "n_cells"))) |>
    dplyr::rename(focal_cell_type = cell_type, focal_type_cells = n_cells)

  pairwise |>
    dplyr::left_join(focal_totals, by = c(strata_cols, "stratum", "focal_cell_type")) |>
    dplyr::left_join(abundance, by = c(strata_cols, "stratum", "neighbor_cell_type")) |>
    dplyr::left_join(focal_abundance, by = c(strata_cols, "stratum", "focal_cell_type")) |>
    dplyr::mutate(
      available_neighbor_cells = stratum_total_cells,
      expected_fraction = ifelse(available_neighbor_cells > 0, neighbor_type_cells / available_neighbor_cells, NA_real_),
      observed_fraction = observed_links / total_links_from_focal_type,
      expected_links = total_links_from_focal_type * expected_fraction,
      enrichment = observed_fraction / expected_fraction
    ) |>
    dplyr::arrange(dplyr::across(dplyr::all_of(c(strata_cols, "focal_cell_type", "neighbor_cell_type"))))
}

tidyr_coalesce_zero <- function(x) {
  x[is.na(x)] <- 0L
  x
}

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    ".cell_type",
    "abundance_fraction",
    "available_neighbor_cells",
    "cell_type",
    "comparison_value",
    "distance",
    "enrichment",
    "expected_fraction",
    "expected_links",
    "reference_value",
    "focal_cell_id",
    "focal_cell_type",
    "focal_type_cells",
    "median_distance",
    "n_cells",
    "neighbor_cell_type",
    "neighbor_type_cells",
    "observed_fraction",
    "observed_links",
    "retained_neighbor_links",
    "stratum",
    "stratum_total_cells",
    "total_links_from_focal_type"
  ))
}
