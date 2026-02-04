#' Estimate cell cycle proliferation index in a monocle3 CDS
#'
#' Computes S phase, G2M phase, and combined proliferation scores using
#' canonical cell cycle gene signatures. Automatically matches both
#' human (upper-case) and mouse (TitleCase) gene symbols.
#'
#' @param cds monocle3 cell_data_set
#'
#' @return CDS with cell cycle scores added to colData
#' @export
estimate_cell_cycle <- function(cds) {

  ## -----------------------------
  ## Canonical cell cycle markers
  ## -----------------------------
  s_genes_human <- c(
    "MCM5","PCNA","TYMS","FEN1","MCM2","MCM4","RRM1","UNG",
    "GINS2","MCM6","CDCA7","DTL","PRIM1","UHRF1","CENPU",
    "HELLS","RFC2","RPA2","NASP","RAD51AP1","GMNN",
    "WDR76","SLBP","CCNE2","UBR7","POLD3","MSH2",
    "ATAD2","RAD51","RRM2","CDC45","CDC6","EXO1",
    "TIPIN","DSCC1","BLM","CASP8AP2","USP1","CLSPN",
    "POLA1","CHAF1B","BRIP1","E2F8",
    # newer additions (keep originals)
    "MCM7","MLF1IP","POLR1B","MRPL36"
  )

  g2m_genes_human <- c(
    "HMGB2","CDK1","NUSAP1","UBE2C","BIRC5","TPX2","TOP2A",
    "NDC80","CKS2","NUF2","CKS1B","MKI67","TMPO","CENPF",
    "TACC3","FAM64A","SMC4","CCNB2","CKAP2L","CKAP2",
    "AURKB","BUB1","KIF11","ANP32E","TUBB4B","GTSE1",
    "KIF20B","HJURP","CDC20","TTK","CDC25C","KIF2C",
    "RANGAP1","NCAPD2","DLGAP5","CDCA3","HMMR",
    "AURKA","PSRC1","ANLN","LBR","CKAP5","CENPE",
    "CTCF","NEK2","G2E3","GAS2L3","CBX5","CENPA","HN1",
    # newer additions
    "PIMREG","JPT1"
  )

  ## --------------------------------
  ## Expand to include mouse symbols
  ## --------------------------------
  to_mouse <- function(x) {
    vapply(x, function(s) {
      s <- tolower(s)
      paste0(toupper(substr(s, 1, 1)), substr(s, 2, nchar(s)))
    }, character(1))
  }

  s_genes <- unique(c(s_genes_human, to_mouse(s_genes_human)))
  g2m_genes <- unique(c(g2m_genes_human, to_mouse(g2m_genes_human)))

  ## --------------------------------
  ## Resolve gene identifiers in CDS
  ## --------------------------------
  if ("gene_short_name" %in% colnames(SummarizedExperiment::rowData(cds))) {
    gene_names <- SummarizedExperiment::rowData(cds)$gene_short_name
  } else {
    gene_names <- rownames(cds)
  }

  gene_names <- as.character(gene_names)

  # counts_mat <- as.matrix(monocle3::counts(cds))
  # bug found: counts not exported by monocle3
  # bug fix: use monocle3::exprs instead
  counts_mat <- Matrix(monocle3::exprs(cds), sparse = TRUE)
  size_factors <- SummarizedExperiment::colData(cds)$Size_Factor

  ## --------------------------------
  ## S phase score
  ## --------------------------------
  s_idx <- gene_names %in% s_genes
  if (any(s_idx)) {
    s_expr <- counts_mat[s_idx, , drop = FALSE]
    # s_expr <- t(t(s_expr) / size_factors)
    # bug found: sparse matrix calculation encountering error
    # bug fix: rewrite calculation for sparse matrix
    s_expr <- s_expr %*%  Diagonal(x = 1 / size_factors)
    s_score <- Matrix::colSums(s_expr)
  } else {
    s_score <- rep(0, ncol(cds))
  }

  ## --------------------------------
  ## G2M phase score
  ## --------------------------------
  g2m_idx <- gene_names %in% g2m_genes
  if (any(g2m_idx)) {
    g2m_expr <- counts_mat[g2m_idx, , drop = FALSE]
    # g2m_expr <- t(t(g2m_expr) / size_factors)
    # bug found: sparse matrix calculation encountering error
    # bug fix: rewrite calculation for sparse matrix
    g2m_expr <- g2m_expr %*% Diagonal(x = 1 / size_factors)
    g2m_score <- Matrix::colSums(g2m_expr)
  } else {
    g2m_score <- rep(0, ncol(cds))
  }

  ## --------------------------------
  ## Store results
  ## --------------------------------
  SummarizedExperiment::colData(cds)$g1s_score <- log1p(s_score)
  SummarizedExperiment::colData(cds)$g2m_score <- log1p(g2m_score)
  SummarizedExperiment::colData(cds)$proliferation_index <-
    log1p(s_score + g2m_score)

  # add a gc() to release memory usage due to dense matrix
  gc()
  return(cds)
}
