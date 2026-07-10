#' Build and preprocess a combined Visium monocle3 cell_data_set
#'
#' Loads multiple Visium samples from CellRanger h5 files,
#' combines them into a single monocle3 cell_data_set, attaches
#' sample metadata, performs basic QC filtering, and runs standard
#' preprocessing and clustering.
#'
#' @param sample_table A data.frame or tibble produced by
#'   \code{discover_visium_sample_table()}, containing at least:
#'   \code{sample_id}, \code{sample_label}, \code{sample_group},
#'   \code{h5_file}, and \code{barcode_suffix}.
#' @param min_umi Minimum total UMIs per spot.
#' @param min_genes Minimum detected genes per spot.
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
build_visium_cds <- function(
    sample_table,
    min_umi     = 100,
    min_genes   = 10,
    num_dim     = 25,
    k           = 20,
    cluster_res = NULL,
    random_seed = 12345
) {

  required_cols <- c(
    "sample_id",
    "sample_label",
    "sample_group",
    "h5_file",
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
  ## Load each Visium h5 into a CDS
  ## ------------------------------------------------------------
  cds_list <- vector("list", length = n)

  for (i in seq_len(n)) {
    message("Reading Visium h5: ", sample_table$h5_file[i])

    cds_i <- read.cds.cellranger.h5.file(sample_table$h5_file[i])

    # Enforce unique barcodes across samples
    colnames(cds_i) <- paste0(colnames(cds_i), sample_table$barcode_suffix[i])

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
  ## Basic QC filtering
  ## ------------------------------------------------------------
  total_umi <- Matrix::colSums(monocle3::exprs(cds))
  cds <- cds[, total_umi > min_umi]

  cds <- monocle3::detect_genes(cds)
  cds <- monocle3::estimate_size_factors(cds)

  cds <- cds[
    ,
    SummarizedExperiment::colData(cds)$num_genes_expressed > min_genes
  ]

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
