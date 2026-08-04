#' Plot an MA expression plot
#'
#' Create an MA plot comparing the mean expression of shared genes between two
#' expression matrices. The most highly expressed genes in each dataset are
#' highlighted and labelled.
#'
#' @param expr1 A matrix-like object containing expression values.
#' @param expr2 A matrix-like object containing expression values.
#' @param shared_genes Character vector of genes present in both datasets.
#' @param top_n Number of highly expressed genes to label from each dataset.
#' @param labels Character vector of length two giving the dataset names.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' plot_ma_expression(expr_mouse, expr_human, shared_genes)
#'
#' @export
plot_ma_expression <- function(
    expr1,
    expr2,
    shared_genes,
    top_n = 50,
    labels = c("Dataset 1", "Dataset 2")
) {

    if (length(labels) != 2) {
        stop("'labels' must have length 2.", call. = FALSE)
    }

    if (!is.numeric(top_n) || length(top_n) != 1 || top_n < 1) {
        stop("'top_n' must be a positive integer.", call. = FALSE)
    }

    shared_genes <- intersect(
        shared_genes,
        intersect(rownames(expr1), rownames(expr2))
    )

    if (length(shared_genes) == 0) {
        stop(
            "No shared genes found in both expression matrices.",
            call. = FALSE
        )
    }

    expr1_mean <- Matrix::rowMeans(
        expr1[shared_genes, , drop = FALSE]
    )

    expr2_mean <- Matrix::rowMeans(
        expr2[shared_genes, , drop = FALSE]
    )

    df <- data.frame(
        Gene = shared_genes,
        Expression1 = expr1_mean,
        Expression2 = expr2_mean
    )

    ## MA transformation
    df$A <- (df$Expression1 + df$Expression2) / 2
    df$M <- df$Expression2 - df$Expression1

    top_expr1 <- head(
        df$Gene[order(df$Expression1, decreasing = TRUE)],
        top_n
    )

    top_expr2 <- head(
        df$Gene[order(df$Expression2, decreasing = TRUE)],
        top_n
    )

    df$Highlight <- df$Gene %in% union(top_expr1, top_expr2)

    ggplot2::ggplot(
        df,
        ggplot2::aes(A, M)
    ) +
        ggplot2::geom_point(
            ggplot2::aes(colour = Highlight),
            alpha = 0.5,
            size = 1.2
        ) +
        ggplot2::scale_color_manual(
            values = c(
                "FALSE" = "grey70",
                "TRUE" = "red"
            ),
            guide = "none"
        ) +
        ggplot2::geom_hline(
            yintercept = 0,
            colour = "black"
        ) +
        ggrepel::geom_text_repel(
            data = df[df$Highlight, ],
            ggplot2::aes(label = Gene),
            size = 3,
            max.overlaps = Inf
        ) +
        ggplot2::theme_classic() +
        ggplot2::labs(
            x = paste("Average expression (", labels[1], " + ", labels[2], ")", sep = ""),
            y = paste(labels[2], "-", labels[1], "expression")
        )
}