#' Plot quality control metrics
#'
#' Generate quality control histograms for single-cell RNA-seq datasets,
#' including UMI counts, mitochondrial percentage, and Scrublet scores.
#' QC thresholds are used to calculate the number of retained cells.
#'
#' @param qc_df A data.frame containing QC metrics.
#' @param n_umi_col Character name of the column containing UMI counts per cell.
#' @param percent_mito_col Character name of the column containing mitochondrial
#'   percentage values.
#' @param scrublet_col Character name of the column containing Scrublet scores.
#' @param umi_min Optional minimum UMI threshold.
#' @param umi_max Optional maximum UMI threshold.
#' @param mito_max Optional maximum mitochondrial percentage threshold.
#' @param scrublet_max Optional maximum Scrublet score threshold.
#' @param bins Number of histogram bins.
#' @param return_list Logical. If TRUE, return plots together with QC filtering
#'   results. If FALSE, return only the plots.
#'
#' @return If `return_list = TRUE`, a named list containing:
#' \describe{
#'   \item{umi}{Histogram of UMI counts.}
#'   \item{mito}{Histogram of mitochondrial percentage.}
#'   \item{scrublet}{Histogram of Scrublet scores.}
#'   \item{keep}{Logical vector indicating cells passing QC thresholds.}
#'   \item{n_cells_before}{Number of cells before filtering.}
#'   \item{n_cells_after}{Number of cells passing QC thresholds.}
#' }
#'
#' If `return_list = FALSE`, a named list containing:
#' \describe{
#'   \item{umi}{Histogram of UMI counts.}
#'   \item{mito}{Histogram of mitochondrial percentage.}
#'   \item{scrublet}{Histogram of Scrublet scores.}
#' }
#'
#' @examples
#' qc_plots <- plot_qc_metrics(
#'     qc_df,
#'     n_umi_col = "nUMI",
#'     percent_mito_col = "percent_mito",
#'     scrublet_col = "scrublet_score",
#'     umi_min = 500,
#'     umi_max = 50000,
#'     mito_max = 10,
#'     scrublet_max = 0.1
#' )
#'
#' qc_results <- plot_qc_metrics(
#'     qc_df,
#'     n_umi_col = "nUMI",
#'     percent_mito_col = "percent_mito",
#'     scrublet_col = "scrublet_score",
#'     return_list = TRUE
#' )
#'
#' @export
plot_qc_metrics <- function(
    qc_df,
    n_umi_col,
    percent_mito_col,
    scrublet_col,
    umi_min = NULL,
    umi_max = NULL,
    mito_max = NULL,
    scrublet_max = NULL,
    bins = 100,
    return_list = FALSE
) {

    ## ------------------------
    ## Input validation
    ## ------------------------

    if (!is.data.frame(qc_df)) {
        stop(
            "'qc_df' must be a data.frame.",
            call. = FALSE
        )
    }

    required_cols <- c(
        n_umi_col,
        percent_mito_col,
        scrublet_col
    )

    missing_cols <- setdiff(
        required_cols,
        colnames(qc_df)
    )

    if (length(missing_cols) > 0) {
        stop(
            "Missing QC columns: ",
            paste(missing_cols, collapse = ", "),
            call. = FALSE
        )
    }

    if (!is.numeric(bins) ||
        length(bins) != 1 ||
        bins < 1) {
        stop(
            "'bins' must be a positive number.",
            call. = FALSE
        )
    }


    ## ------------------------
    ## Calculate QC filtering
    ## ------------------------

    n_cells_before <- nrow(qc_df)

    keep <- rep(TRUE, n_cells_before)

    if (!is.null(umi_min)) {
        keep <- keep &
            qc_df[[n_umi_col]] >= umi_min
    }

    if (!is.null(umi_max)) {
        keep <- keep &
            qc_df[[n_umi_col]] <= umi_max
    }

    if (!is.null(mito_max)) {
        keep <- keep &
            qc_df[[percent_mito_col]] <= mito_max
    }

    if (!is.null(scrublet_max)) {
        keep <- keep &
            qc_df[[scrublet_col]] <= scrublet_max
    }

    n_cells_after <- sum(
        keep,
        na.rm = TRUE
    )


    cell_text <- paste0(
        "Cells before QC: ",
        n_cells_before,
        "\n",
        "Cells after QC:  ",
        n_cells_after
    )


    ## ------------------------
    ## UMI plot
    ## ------------------------

    p_umi <- ggplot2::ggplot(
        qc_df,
        ggplot2::aes(
            x = .data[[n_umi_col]]
        )
    ) +
        ggplot2::geom_histogram(
            bins = bins,
            fill = "steelblue",
            color = "black"
        ) +
        ggplot2::theme_classic() +
        ggplot2::labs(
            title = "UMI counts per cell",
            x = n_umi_col,
            y = "Cell count"
        ) +
        ggplot2::annotate(
            "text",
            x = Inf,
            y = Inf,
            label = cell_text,
            hjust = 1.1,
            vjust = 1.1,
            size = 4
        )

    if (!is.null(umi_min)) {
        p_umi <- p_umi +
            ggplot2::geom_vline(
                xintercept = umi_min,
                color = "red",
                linetype = "dashed"
            )
    }

    if (!is.null(umi_max)) {
        p_umi <- p_umi +
            ggplot2::geom_vline(
                xintercept = umi_max,
                color = "red",
                linetype = "dashed"
            )
    }


    ## ------------------------
    ## Mitochondrial plot
    ## ------------------------

    p_mito <- ggplot2::ggplot(
        qc_df,
        ggplot2::aes(
            x = .data[[percent_mito_col]]
        )
    ) +
        ggplot2::geom_histogram(
            bins = bins,
            fill = "darkorange",
            color = "black"
        ) +
        ggplot2::theme_classic() +
        ggplot2::labs(
            title = "Mitochondrial UMIs (%)",
            x = percent_mito_col,
            y = "Cell count"
        )

    if (!is.null(mito_max)) {
        p_mito <- p_mito +
            ggplot2::geom_vline(
                xintercept = mito_max,
                color = "red",
                linetype = "dashed"
            )
    }


    ## ------------------------
    ## Scrublet plot
    ## ------------------------

    p_scrub <- ggplot2::ggplot(
        qc_df,
        ggplot2::aes(
            x = .data[[scrublet_col]]
        )
    ) +
        ggplot2::geom_histogram(
            bins = bins,
            fill = "purple",
            color = "black"
        ) +
        ggplot2::theme_classic() +
        ggplot2::labs(
            title = "Scrublet scores",
            x = scrublet_col,
            y = "Cell count"
        )

    if (!is.null(scrublet_max)) {
        p_scrub <- p_scrub +
            ggplot2::geom_vline(
                xintercept = scrublet_max,
                color = "red",
                linetype = "dashed"
            )
    }


    ## ------------------------
    ## Return
    ## ------------------------

    if (return_list) {

        return(
            list(
                umi = p_umi,
                mito = p_mito,
                scrublet = p_scrub,
                keep = keep,
                n_cells_before = n_cells_before,
                n_cells_after = n_cells_after
            )
        )

    } else {

        return(
            list(
                umi = p_umi,
                mito = p_mito,
                scrublet = p_scrub
            )
        )
    }
}