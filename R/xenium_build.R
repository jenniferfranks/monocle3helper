#' Build and preprocess a combined Xenium monocle3 cell_data_set
#'
#' Loads multiple Xenium samples from Xenium Onboard Analysis h5 files,
#' combines them into a single monocle3 cell_data_set, attaches sample
#' and cell-centroid metadata, performs basic QC filtering, and runs
#' standard preprocessing and clustering.
#'
#' QC defaults (minimum total counts, minimum detected genes, and
#' nucleus area bounds) follow Vannan et al.
#'
#' @param sample_table A data.frame or tibble produced by
#'   \code{discover_xenium_sample_table()}, containing at least:
#'   \code{sample_id}, \code{sample_label}, \code{sample_group},
#'   \code{h5_file}, \code{cells_csv}, and \code{barcode_suffix}.
#' @param min_counts Minimum total transcript counts per cell.
#' @param min_genes Minimum detected genes per cell.
#' @param min_nucleus_area Minimum nucleus area (square microns) per cell.
#' @param max_nucleus_area Maximum nucleus area (square microns) per cell.
#' @param num_dim Number of dimensions for PCA preprocessing.
#' @param k Number of nearest neighbors used to build the graph for
#'   Leiden clustering. Larger k gives fewer, coarser clusters.
#' @param cluster_res Optional resolution parameter for clustering. Defaults
#'   to NULL, letting monocle3 determine resolution automatically from k;
#'   set this explicitly only if you need to override that.
#' @param random_seed Random seed for reproducibility.
#'
#' @return A processed monocle3 cell_data_set containing all samples.
#'
#' @export
build_xenium_cds <- function(
    sample_table,
    min_counts       = 12,
    min_genes        = 10,
    min_nucleus_area = 6,
    max_nucleus_area = 80,
    num_dim          = 50,
    k                = 20,
    cluster_res      = NULL,
    random_seed      = 12345
) {

  required_cols <- c(
    "sample_id",
    "sample_label",
    "sample_group",
    "h5_file",
    "cells_csv",
    "barcode_suffix"
  )

  if (!all(required_cols %in% colnames(sample_table))) {
    stop(
      "sample_table must contain columns: ",
      paste(required_cols, collapse = ", ")
    )
  }

  n <- nrow(sample_table)
  if (n == 0) {
    stop("sample_table contains zero rows.")
  }

  ## ------------------------------------------------------------
  ## Load each Xenium h5 into a CDS
  ## ------------------------------------------------------------
  cds_list <- vector("list", length = n)

  for (i in seq_len(n)) {
    message("Reading Xenium h5: ", sample_table$h5_file[i])

    cds_i <- read.cds.cellranger.h5.file(sample_table$h5_file[i])

    # Track the original h5 barcode and originating sample via colData:
    # monocle3::combine_cds() applies its own disambiguating suffix to
    # colnames() when merging, so we can't rely on colnames() alone to
    # reconstruct a stable cell ID after combining.
    SummarizedExperiment::colData(cds_i)$raw_barcode <- colnames(cds_i)
    SummarizedExperiment::colData(cds_i)$sample_numeric <- i

    cds_list[[i]] <- cds_i
  }

  ## ------------------------------------------------------------
  ## Combine CDS objects
  ## ------------------------------------------------------------
  cds <- monocle3::combine_cds(cds_list)

  ## ------------------------------------------------------------
  ## Map sample metadata
  ## ------------------------------------------------------------
  numeric_sample <- as.character(
    SummarizedExperiment::colData(cds)$sample_numeric
  )

  # Rebuild colnames from the tracked raw barcode + our own barcode_suffix,
  # overriding whatever combine_cds() renamed them to, so downstream joins
  # against cells_csv (which use this same barcode_suffix convention) work.
  colnames(cds) <- paste0(
    SummarizedExperiment::colData(cds)$raw_barcode,
    sample_table$barcode_suffix[as.integer(numeric_sample)]
  )

  map_sample_id <- stats::setNames(
    sample_table$sample_id,
    seq_len(n)
  )

  map_sample_label <- stats::setNames(
    sample_table$sample_label,
    seq_len(n)
  )

  map_sample_group <- stats::setNames(
    sample_table$sample_group,
    seq_len(n)
  )

  SummarizedExperiment::colData(cds)$sample <- unname(
    map_sample_id[numeric_sample]
  )

  SummarizedExperiment::colData(cds)$sample_label <- unname(
    map_sample_label[numeric_sample]
  )

  SummarizedExperiment::colData(cds)$sample_group <- unname(
    map_sample_group[numeric_sample]
  )

  ## ------------------------------------------------------------
  ## Remove deprecated / control probes
  ## ------------------------------------------------------------
  cds <- cds[!grepl("Deprecated", rownames(cds)), ]
  cds <- cds[!grepl("Codeword", rownames(cds)), ]
  cds <- cds[!grepl("NegControl", rownames(cds)), ]

  ## ------------------------------------------------------------
  ## Basic QC filtering
  ## ------------------------------------------------------------
  total_counts <- Matrix::colSums(monocle3::exprs(cds))
  cds <- cds[, total_counts > min_counts]

  cds <- monocle3::detect_genes(cds)
  cds <- monocle3::estimate_size_factors(cds)

  cds <- cds[
    ,
    SummarizedExperiment::colData(cds)$num_genes_expressed > min_genes
  ]

  ## ------------------------------------------------------------
  ## Attach cell centroid metadata
  ## ------------------------------------------------------------
  cell_meta_list <- vector("list", length = n)

  for (i in seq_len(n)) {
    cell_meta_list[[i]] <- read_xenium_cells(
      cells_csv      = sample_table$cells_csv[i],
      barcode_suffix = sample_table$barcode_suffix[i]
    )
  }

  cell_meta <- dplyr::bind_rows(cell_meta_list)

  cd <- as.data.frame(SummarizedExperiment::colData(cds))
  cd$cell_id <- rownames(cd)

  result <- dplyr::left_join(cd["cell_id"], cell_meta, by = "cell_id")

  if (any(rownames(SummarizedExperiment::colData(cds)) != result$cell_id)) {
    stop("Cell ID mismatch while attaching Xenium spatial metadata.")
  }

  meta_cols <- setdiff(colnames(result), "cell_id")
  for (col in meta_cols) {
    SummarizedExperiment::colData(cds)[[col]] <- result[[col]]
  }

  n_unmatched <- sum(is.na(result$x_centroid))
  if (n_unmatched > 0) {
    warning(
      n_unmatched, " of ", nrow(result),
      " cells had no matching row in cells_csv and will be dropped ",
      "(check that cells_csv corresponds to the same Xenium run as h5_file)."
    )
    cds <- cds[, !is.na(SummarizedExperiment::colData(cds)$x_centroid)]
  }

  ## ------------------------------------------------------------
  ## QC filtering on nucleus area
  ## ------------------------------------------------------------
  if ("nucleus_area" %in% colnames(SummarizedExperiment::colData(cds))) {
    nucleus_area <- SummarizedExperiment::colData(cds)$nucleus_area
    cds <- cds[
      ,
      !is.na(nucleus_area) &
        nucleus_area > min_nucleus_area &
        nucleus_area < max_nucleus_area
    ]
  }

  ## ------------------------------------------------------------
  ## Preprocessing and clustering
  ## ------------------------------------------------------------
  set.seed(random_seed)

  cds <- monocle3::preprocess_cds(
    cds,
    num_dim = num_dim
  )

  cds <- monocle3::align_cds(
    cds,
    alignment_group = "sample"
  )

  cds <- monocle3::reduce_dimension(
    cds,
    umap.fast_sgd = TRUE
  )

  cds <- monocle3::cluster_cells(
    cds,
    k = k,
    resolution = cluster_res,
    random_seed = random_seed
  )

  SummarizedExperiment::colData(cds)$clusters <-
    monocle3::clusters(cds)

  cds
}
