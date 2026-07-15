test_that("xenium proximity reports directed pairwise enrichment", {
  skip_if_not_installed("RANN")
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("S4Vectors")

  counts <- matrix(0, nrow = 1, ncol = 5)
  colnames(counts) <- paste0("cell", seq_len(5))
  cds <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    colData = S4Vectors::DataFrame(
      cell_type = c("A", "B", "B", "A", "B"),
      sample = c("s1", "s1", "s1", "s2", "s2"),
      roi_id = c("r1", "r1", "r1", "r2", "r2"),
      x_centroid = c(0, 1, 10, 0, 2),
      y_centroid = c(0, 0, 0, 0, 0)
    )
  )

  proximity <- analyze_xenium_cell_proximity(
    cds,
    cell_type_col = "cell_type",
    strata_cols = "sample",
    k = 1,
    max_distance = 5
  )

  expect_s3_class(proximity, "xenium_cell_proximity")
  expect_true(all(c("pairwise", "by_stratum", "cell_neighbors") %in% names(proximity)))
  expect_true(all(c("focal_cell_type", "neighbor_cell_type", "enrichment", "median_distance") %in%
    colnames(proximity$pairwise)))
  expect_true(all(proximity$cell_neighbors$stratum %in% c("s1", "s2")))
})

test_that("xenium proximity can stratify by ROI and exclude same-type links", {
  skip_if_not_installed("RANN")
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("S4Vectors")

  counts <- matrix(0, nrow = 1, ncol = 4)
  colnames(counts) <- paste0("cell", seq_len(4))
  cds <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    colData = S4Vectors::DataFrame(
      cell_type = c("A", "A", "B", "B"),
      sample = "s1",
      roi_id = c("r1", "r1", "r1", "r2"),
      x_centroid = c(0, 1, 2, 100),
      y_centroid = c(0, 0, 0, 0)
    )
  )

  proximity <- analyze_xenium_cell_proximity(
    cds,
    cell_type_col = "cell_type",
    strata_cols = "sample",
    roi_col = "roi_id",
    k = 1,
    max_distance = Inf,
    exclude_same_type = TRUE
  )

  expect_true("roi_id" %in% colnames(proximity$pairwise))
  expect_false(any(proximity$cell_neighbors$focal_cell_type == proximity$cell_neighbors$neighbor_cell_type))
})


test_that("xenium proximity compares sample-level statistics", {
  skip_if_not_installed("RANN")
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("S4Vectors")

  counts <- matrix(0, nrow = 1, ncol = 6)
  colnames(counts) <- paste0("cell", seq_len(6))
  cds <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    colData = S4Vectors::DataFrame(
      cell_type = c("A", "B", "B", "A", "B", "A"),
      sample = c("control", "control", "control", "treated", "treated", "treated"),
      roi_id = c("r1", "r1", "r1", "r1", "r1", "r1"),
      x_centroid = c(0, 1, 10, 0, 5, 6),
      y_centroid = c(0, 0, 0, 0, 0, 0)
    )
  )

  proximity <- analyze_xenium_cell_proximity(
    cds,
    cell_type_col = "cell_type",
    strata_cols = c("sample", "roi_id"),
    k = 1,
    max_distance = Inf
  )

  differences <- compare_xenium_cell_proximity(
    proximity,
    stratum_col = "sample",
    within_cols = "roi_id",
    reference = "control",
    value = "median_distance"
  )

  expect_true(all(c("sample", "reference", "roi_id", "difference", "log2_ratio") %in% colnames(differences)))
  expect_true(all(differences$sample == "treated"))
})
