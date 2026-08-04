#' Plot mean expression scatter
#'
#' Compare mean expression between two expression matrices for a shared set
#' of genes. The most highly expressed genes in each dataset are highlighted
#' and labelled.
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
#' plot_mean_expression_scatter(expr_mouse, expr_human, shared_genes)
#'
#' @export
plot_mean_expression_scatter <- function(
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

    rho <- stats::cor(
        expr1_mean,
        expr2_mean,
        method = "spearman",
        use = "complete.obs"
    )

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
        ggplot2::aes(Expression1, Expression2)
    ) +
        ggplot2::geom_point(
            ggplot2::aes(colour = Highlight),
            alpha = 0.6,
            size = 1.5
        ) +
        ggplot2::scale_color_manual(
            values = c(
                "FALSE" = "grey70",
                "TRUE" = "red"
            ),
            guide = "none"
        ) +
        ggplot2::geom_smooth(
            method = "lm",
            colour = "black",
            se = FALSE
        ) +
        ggrepel::geom_text_repel(
            data = df[df$Highlight, ],
            ggplot2::aes(label = Gene),
            size = 3,
            max.overlaps = Inf
        ) +
        ggplot2::annotate(
            "text",
            x = Inf,
            y = Inf,
            hjust = 1.1,
            vjust = 1.5,
            label = paste0(
                "Spearman = ",
                round(rho, 3)
            )
        ) +
        ggplot2::theme_classic() +
        ggplot2::labs(
            x = paste("Mean", labels[1], "expression"),
            y = paste("Mean", labels[2], "expression")
        )
}