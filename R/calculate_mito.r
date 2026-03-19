#' @export
# Mitochondrial reads filtering 
#' @export
calculate_mito <- function(cds, pattern = "^MT-|^mt-",
                           col_name = "perc_mitochondrial_umis") {
  
  all_genes <- rowData(cds)$gene_short_name
  
  if (is.null(all_genes)) {
    stop("rowData(cds)$gene_short_name is missing.")
  }
  
  all_genes <- as.character(all_genes)
  
  # find mito genes by pattern
  mito_idx <- grepl(pattern, all_genes)
  
  # fallback for alternative naming
  if (!any(mito_idx, na.rm = TRUE)) {
    mito_genes_base <- c(
      "Atp6", "Atp8",
      "Co1", "Co2", "Co3",
      "Cytb",
      "Nd1", "Nd2", "Nd3", "Nd4", "Nd4l", "Nd5", "Nd6"
    )
    
    mito_names <- c(mito_genes_base, toupper(mito_genes_base))
    mito_idx <- all_genes %in% mito_names
  }
  
  # no mito genes found
  if (!any(mito_idx, na.rm = TRUE)) {
    colData(cds)[[col_name]] <- NA_real_
    return(cds)
  }
  
  counts_matrix <- counts(cds)
  
  total_mito_counts <- Matrix::colSums(counts_matrix[mito_idx, , drop = FALSE])
  total_counts <- Matrix::colSums(counts_matrix)
  
  percent_mito <- ifelse(
    total_counts > 0,
    (total_mito_counts / total_counts) * 100,
    0
  )
  
  colData(cds)[[col_name]] <- percent_mito
  
  return(cds)
}
