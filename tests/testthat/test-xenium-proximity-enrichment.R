make_enrichment_test_cds <- function() {
  x_centroid <- c(0, 1, 2, 3, 100, 101, 102, 103)
  y_centroid <- rep(0, 8)
  cell_type <- c(rep("A", 4), rep("B", 4))
  sample <- rep("s1", 8)
  roi_id <- c(rep("r1", 4), rep("r2", 4))

  counts <- matrix(0, nrow = 1, ncol = length(x_centroid))
  colnames(counts) <- paste0("cell", seq_along(x_centroid))
  SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    colData = S4Vectors::DataFrame(
      cell_type = cell_type,
      sample = sample,
      roi_id = roi_id,
      x_centroid = x_centroid,
      y_centroid = y_centroid
    )
  )
}

test_that("xenium_proximity_enrichment reports permutation-based enrichment per sample", {
  skip_if_not_installed("RANN")
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("S4Vectors")

  cds <- make_enrichment_test_cds()

  result <- xenium_proximity_enrichment(
    cds,
    cell_type_col = "cell_type",
    strata_cols = "sample",
    network_method = "knn",
    k = 2,
    n_permutations = 99,
    seed = 42,
    verbose = FALSE
  )

  expect_s3_class(result, "xenium_proximity_enrichment")
  expect_true(all(c("enrichment", "by_stratum") %in% names(result)))
  expect_true(all(c(
    "cell_type_a", "cell_type_b", "observed", "expected", "sd_sim",
    "enrichment", "p_value", "p_adj", "interaction_type"
  ) %in% colnames(result$enrichment)))

  # Two far-apart, purely same-type clusters with a KNN network small enough
  # to never bridge them: cross-type edges must be exactly zero, and both
  # homotypic (same-type) combinations must be present.
  ab_row <- result$enrichment[
    result$enrichment$cell_type_a == "A" & result$enrichment$cell_type_b == "B",
  ]
  expect_equal(nrow(ab_row), 1)
  expect_equal(ab_row$observed, 0)

  expect_true(all(c("A--A", "B--B") %in%
    paste(result$enrichment$cell_type_a, result$enrichment$cell_type_b, sep = "--")))
  expect_true(all(result$enrichment$p_value >= 0 & result$enrichment$p_value <= 1))
  expect_true(all(is.finite(result$enrichment$enrichment)))
  expect_equal(sort(unique(result$by_stratum$cell_type)), c("A", "B"))
})

test_that("xenium_proximity_enrichment can stratify by ROI", {
  skip_if_not_installed("RANN")
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("S4Vectors")

  cds <- make_enrichment_test_cds()

  result <- xenium_proximity_enrichment(
    cds,
    cell_type_col = "cell_type",
    strata_cols = "sample",
    roi_col = "roi_id",
    network_method = "knn",
    k = 2,
    n_permutations = 49,
    seed = 1,
    verbose = FALSE
  )

  expect_true("roi_id" %in% colnames(result$enrichment))
  expect_true("roi_id" %in% colnames(result$by_stratum))
  expect_equal(sort(unique(result$enrichment$roi_id)), c("r1", "r2"))
  # Each ROI here contains only one cell type, so no cross-type pairs exist.
  expect_true(all(result$enrichment$cell_type_a == result$enrichment$cell_type_b))
})

test_that("xenium_proximity_enrichment validates arguments", {
  skip_if_not_installed("RANN")
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("S4Vectors")

  cds <- make_enrichment_test_cds()

  expect_error(
    xenium_proximity_enrichment(cds, network_method = "knn", k = 0),
    "k must be a positive integer"
  )
  expect_error(
    xenium_proximity_enrichment(cds, network_method = "knn", n_permutations = 0),
    "n_permutations must be a positive integer"
  )
  expect_error(
    xenium_proximity_enrichment(cds, network_method = "knn", sig_threshold = 1.5),
    "sig_threshold must be a number between 0 and 1"
  )
})

test_that("Xenium proximity enrichment plots return ggplot objects", {
  skip_if_not_installed("RANN")
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("S4Vectors")
  skip_if_not_installed("ggplot2")

  cds <- make_enrichment_test_cds()
  result <- xenium_proximity_enrichment(
    cds,
    cell_type_col = "cell_type",
    strata_cols = "sample",
    network_method = "knn",
    k = 2,
    n_permutations = 49,
    seed = 7,
    verbose = FALSE
  )

  expect_s3_class(plot_xenium_proximity_heatmap(result), "gg")
  expect_s3_class(
    plot_xenium_proximity_barplot(result, significant_only = FALSE),
    "gg"
  )
  expect_s3_class(
    plot_xenium_proximity_network(result, significant_only = FALSE),
    "gg"
  )
})
