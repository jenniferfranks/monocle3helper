#' @export
# Mitochondrial reads filtering 
calculate_mito <- function(cds, pattern = "^MT-|^mt-") {
  
  # Identify mitochondrial genes based on the pattern
  all_genes <- fData(cds)$gene_short_name
  mito_genes <- grep(pattern = pattern, x = all_genes, value = TRUE)
  

  counts_matrix <- counts(cds)
  
  #  Calculate mitochondrial counts per cell
  mito_counts_matrix <- counts_matrix[mito_genes, , drop = FALSE]
  total_mito_counts <- Matrix::colSums(mito_counts_matrix)
  
  total_counts <- Matrix::colSums(counts_matrix)
  
  # Calculate percentage, handling division by zero (cells with 0 total counts)
  percent_mito <- ifelse(
    total_counts > 0,
    (total_mito_counts / total_counts) * 100,
    0
  ) 
  
  #Add the new column to the cell metadata (colData, which is pData in older SingleCellExperiment/Monocle versions)
  colData(cds)[["perc_mitochondrial_umis"]] <- percent_mito
  
  return(cds)
}
