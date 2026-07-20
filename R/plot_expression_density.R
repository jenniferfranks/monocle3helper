#' Plot expression density
#'
#' Compare the distribution of expression values between two expression
#' matrices. For a single gene, the density is computed across cells. For
#' multiple genes, the density is computed from the mean expression of each
#' gene across all cells.
#'
#' @param expr1 A matrix-like object containing expression values.
#' @param expr2 A matrix-like object containing expression values.
#' @param genes Character vector of gene names.
#' @param labels Character vector of length two giving the labels for
#'   \code{expr1} and \code{expr2}.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' plot_expression_density(expr_mouse, expr_human, "Gad1")
#' plot_expression_density(expr_mouse, expr_human, marker_genes)
#'
#' @export
plot_expression_density <- function(
    expr1,
    expr2,
    genes,
    labels = c("Dataset 1", "Dataset 2")
) {

    if (length(labels) != 2) {
        stop("'labels' must have length 2.", call. = FALSE)
    }

    genes <- intersect(
        genes,
        intersect(rownames(expr1), rownames(expr2))
    )

    if (length(genes) == 0) {
        stop(
            "None of the supplied genes were found in both expression matrices.",
            call. = FALSE
        )
    }

    if (length(genes) == 1) {

        x1 <- as.numeric(expr1[genes, , drop = FALSE])
        x2 <- as.numeric(expr2[genes, , drop = FALSE])

        xlab <- paste0(genes, " expression")

    } else {

        x1 <- Matrix::rowMeans(
            expr1[genes, , drop = FALSE]
        )

        x2 <- Matrix::rowMeans(
            expr2[genes, , drop = FALSE]
        )

        xlab <- "Mean expression"
    }

    df <- rbind(
        data.frame(
            Expression = x1,
            Dataset = labels[1]
        ),
        data.frame(
            Expression = x2,
            Dataset = labels[2]
        )
    )

    ggplot2::ggplot(
        df,
        ggplot2::aes(
            x = Expression,
            colour = Dataset,
            fill = Dataset
        )
    ) +
        ggplot2::geom_density(alpha = 0.3) +
        ggplot2::theme_classic() +
        ggplot2::labs(
            x = xlab,
            y = "Density",
            colour = NULL,
            fill = NULL
        )
}