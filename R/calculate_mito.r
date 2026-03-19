#' @export
# Mitochondrial reads filtering 
#' @export
calculate_mito <- function(cds, pattern = "^MT-|^mt-") {

  # Identify mitochondrial genes based on the pattern
  all_genes <- rowData(cds)$gene_short_name
  mito_genes <- grep(pattern = pattern, x = all_genes, value = TRUE)


  # Handle case where no mitochondrial genes are labeled differently
  if (length(mito_genes) == 0) {
    mito_genes_base <- c(
      "Atp6", "Atp8",
      "Co1", "Co2", "Co3",
      "Cytb",
      "Nd1", "Nd2", "Nd3", "Nd4", "Nd4l", "Nd5", "Nd6"
    )
    mito_genes <- c(mito_genes_base, toupper(mito_genes_base))

    #check if these exist in CDS
    if(!any(mito_genes %in% all_genes)){
      # If no mitochondrial genes are found, assign NA/0 to the new column
      colData(cds)[[col_name]] <- NA
      return(cds)
    }
  }
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
