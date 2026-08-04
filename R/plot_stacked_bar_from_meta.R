#' Plot stacked bar chart from metadata
#' 
plot_stacked_bar_from_meta <- function(
  meta1,             # first metadata column (e.g., cell type)
  meta2,             # second metadata column (e.g., sample or condition)
  col_order = NULL,
  fill_by = "meta1",
  fill_colors = NULL,
  fill_order = NULL,
  fill_title = "legend_title",
  xlab = "cell type",
  ylab = "percentage of cells",
  normalize = c("none", "sample", "celltype", "both")
) {

  normalize <- match.arg(normalize)

  # create contingency table
  tab <- table(meta1, meta2)

  # normalize
  if (normalize == "sample") {
    # proportions within each sample (columns sum to 1)
    tab <- prop.table(tab, margin = 2)

  } else if (normalize == "celltype") {
    # proportions within each cell type (rows sum to 1)
    tab <- prop.table(tab, margin = 1)

  } else if (normalize == "both") {
    # Step 1: normalize by total cells in each sample
    tab <- prop.table(tab, margin = 1)

    # Step 2: normalize each cell type across samples
    tab <- prop.table(tab, margin = 2)
  }

  # reorder columns if requested
  if (!is.null(col_order)) {
    tab <- tab[, col_order, drop = FALSE]
  }

  # convert to long format
  df <- data.frame(
    meta1 = rep(rownames(tab), ncol(tab)),
    meta2 = rep(colnames(tab), each = nrow(tab)),
    value = as.vector(tab)
  )

  # determine x and fill variables
  fill_var <- if (fill_by == "meta1") "meta1" else "meta2"
  x_var    <- if (fill_by == "meta1") "meta2" else "meta1"

  # ordering
  if (!is.null(fill_order)) {
    df[[fill_var]] <- factor(df[[fill_var]], levels = fill_order)
  }

  if (!is.null(col_order)) {
    df[[x_var]] <- factor(df[[x_var]], levels = col_order)
  }

  # plot
  p <- ggplot(df, aes(
    x = .data[[x_var]],
    y = value,
    fill = .data[[fill_var]]
  )) +
    geom_col() +
    {if (!is.null(fill_colors)) scale_fill_manual(values = fill_colors) else NULL} +
    labs(fill = fill_title, x = xlab, y = ylab) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  return(p)
}