#' Score Multiple Gene Sets in a Monocle3 cell_data_set
#'
#' Computes aggregate expression scores for multiple gene sets using
#' \link[monocle3]{aggregate_gene_expression} and returns a cell-level
#' data frame containing the original phenotype data and gene set scores.
#'
#' Each gene set is converted into a binary membership vector indicating
#' whether each gene in the dataset belongs to the set. Expression values
#' are normalized, optionally scaled, and aggregated according to user-
#' specified parameters.
#'
#' @param cds A \link[monocle3]{cell_data_set} object.
#'
#' @param gene_sets A named list of character vectors. Each element
#'   corresponds to a gene set and must contain gene symbols matching
#'   \code{rowData(cds)$gene_short_name}.
#'
#' @param scale_agg Logical. If TRUE, aggregated values are Z-score scaled
#'   across cells.
#'
#' @param norm_method Character string specifying the normalization
#'   method. Passed to \code{aggregate_gene_expression()}.
#'
#' @param gene_agg_fun Character string specifying how gene-level expression
#'   is aggregated. Must be one of \code{"sum"} or \code{"mean"}.
#'   Default is \code{"sum"}.
#'
#' @return A data.frame containing \code{pData(cds)} with additional
#'   columns corresponding to gene set scores.
#'
#'#' @examples
#' \dontrun{
#' gene_sets <- list(
#'   TGFb = c("Tgfb1", "Smad3", "Smad2"),
#'   EMT  = c("Vim", "Cdh2", "Zeb1")
#' )
#'
#'
#'
#' scores <- score_gene_sets(
#'   cds,
#'   gene_sets,
#'   scale_agg = TRUE,
#'   norm_method = "log",
#'   gene_agg_fun = "mean"
#' )
#' }
#'
#' @export
#'
#'
score_gene_sets <- function(cds,
                            gene_sets,
                            scale_agg = TRUE,
                            norm_method = "log",
                            gene_agg_fun = c("sum", "mean")) {

  # Match aggregation function
  gene_agg_fun <- match.arg(gene_agg_fun)

  # Ensure gene_sets is a named list
  if (!is.list(gene_sets) || is.null(names(gene_sets))) {
    stop("gene_sets must be a named list of character vectors.")
  }

  # Extract phenotype data
  results <- pData(cds)

  # Extract gene symbols
  gene_names <- rowData(cds)$gene_short_name

  if (is.null(gene_names)) {
    stop("rowData(cds)$gene_short_name is missing.")
  }

  gene_names <- as.character(gene_names)

  # Loop over gene sets
  for (i in seq_along(gene_sets)) {

    set_name <- names(gene_sets)[i]
    gene_list <- unique(gene_sets[[i]])


    gene.group <- data.frame(
      gene  = gene_names,
      group = ifelse(gene_names %in% gene_list, set_name, "other"),
      row.names = rownames(rowData(cds))
    )

    agg <- aggregate_gene_expression(
      cds, gene.group,
      scale_agg_values = scale_agg,
      norm_method = norm_method,
      gene_agg_fun = gene_agg_fun
    )

    if (!(set_name %in% rownames(agg))) next
    results[[set_name]] <- as.numeric(agg[set_name, ])
  }

  return(results)
}
