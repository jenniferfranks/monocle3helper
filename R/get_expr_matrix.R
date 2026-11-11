#' Retrieve an expression matrix
#'
#' Extract an expression matrix from a
#' \code{SingleCellExperiment}-derived object (including Monocle3
#' \code{cell_data_set} objects). If the requested assay is
#' `"logcounts"` but no log-normalized assay is present, a
#' log-transformed count matrix is generated on the fly using
#' `log1p(counts(cds))`.
#'
#' @param cds A \code{SingleCellExperiment} or Monocle3
#'   \code{cell_data_set} object.
#' @param assay Name of the assay to retrieve.
#'
#' @return A matrix-like object containing expression values.
#'
#' @examples
#' expr <- get_expr_matrix(cds)
#' expr <- get_expr_matrix(cds, assay = "counts")
#'
#' @export
get_expr_matrix <- function(
    cds,
    assay = "logcounts"
) {

    if (!inherits(cds, "SingleCellExperiment")) {
        stop(
            "'cds' must inherit from 'SingleCellExperiment'.",
            call. = FALSE
        )
    }

    available_assays <- SummarizedExperiment::assayNames(cds)

    if (assay == "logcounts") {

        if ("logcounts" %in% available_assays) {
            return(
                SummarizedExperiment::assay(
                    cds,
                    "logcounts"
                )
            )
        }

        warning(
            "No 'logcounts' assay found. Returning log1p(counts(cds)).",
            call. = FALSE
        )

        return(
            log1p(
                SingleCellExperiment::counts(cds)
            )
        )
    }

    if (!assay %in% available_assays) {
        stop(
            "Assay '", assay,
            "' not found. Available assays: ",
            paste(available_assays, collapse = ", "),
            call. = FALSE
        )
    }

    SummarizedExperiment::assay(cds, assay)
}