#' Run Scrublet doublet detection on a CDS
#'
#' @param cds monocle3 cell_data_set
#' @param expected_doublet_rate Numeric
#' @param conda_env Conda environment name
#'
#' @return CDS with scrublet_score and scrublet_call
#' @export
run_scrublet <- function(
    cds,
    expected_doublet_rate = 0.06,
    conda_env = NULL
) {

  if (!is.null(conda_env)) {
    reticulate::use_condaenv(conda_env, required = TRUE)
  }

  counts_matrix <- t(monocle3::exprs(cds))

  reticulate::py_run_string("import scrublet as scr")

  scrub <- reticulate::py$scr$Scrublet(
    counts_matrix,
    expected_doublet_rate = expected_doublet_rate
  )

  res <- scrub$scrub_doublets(
    min_counts = 2,
    min_cells = 3,
    min_gene_variability_pctl = 85,
    n_prin_comps = 30
  )

  SummarizedExperiment::pData(cds)$scrublet_score <- as.numeric(res[[1]])
  SummarizedExperiment::pData(cds)$scrublet_call  <- res[[2]]

  cds
}
