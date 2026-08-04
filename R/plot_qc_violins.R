#' Plot QC metrics as violin plots
#' 
plot_qc_violins <- function(
    metadata_df,
    metrics,
    metric_labels = NULL,
    plot_titles = NULL,
    group_col = NULL,
    y_limits = NULL,
    violin_color = "#4C72B0",
    point_size = 0,
    base_size = 18
) {

  library(ggplot2)
  library(dplyr)
  library(tidyr)

  ## Check columns exist
  missing_cols <- setdiff(metrics, colnames(metadata_df))
  if (length(missing_cols) > 0) {
    stop("Missing columns: ", paste(missing_cols, collapse = ", "))
  }

  ## Defaults
  if (is.null(metric_labels))
    metric_labels <- metrics

  if (is.null(plot_titles))
    plot_titles <- metric_labels

  if (length(metric_labels) != length(metrics))
    stop("metric_labels must have the same length as metrics.")

  if (length(plot_titles) != length(metrics))
    stop("plot_titles must have the same length as metrics.")

  ## Long format
  if (!is.null(group_col)) {

    plot_df <- metadata_df %>%
      select(all_of(c(group_col, metrics))) %>%
      pivot_longer(
        cols = all_of(metrics),
        names_to = "Metric",
        values_to = "Value"
      )

    x_var <- group_col

  } else {

    plot_df <- metadata_df %>%
      select(all_of(metrics)) %>%
      pivot_longer(
        everything(),
        names_to = "Metric",
        values_to = "Value"
      )

    plot_df$Group <- ""
    x_var <- "Group"
  }

  ## Rename facets
  plot_df$Metric <- factor(
    plot_df$Metric,
    levels = metrics,
    labels = plot_titles
  )

  ## Named vector for x-axis labels
  #labeller <- setNames(metric_labels, plot_titles)

  p <- ggplot(
    plot_df,
    aes(
      x = .data[[x_var]],
      y = Value
    )
  ) +
    geom_violin(
      fill = violin_color,
      color = "black",
      trim = FALSE,
      linewidth = 0.6
    ) +
    facet_wrap(
      ~Metric,
      scales = "free_y",
      #labeller = labeller(Metric = labeller)
    ) +
    theme_classic(base_size = base_size) +
    labs(
      x = NULL,
      y = NULL
    ) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(
        face = "bold",
        size = base_size + 2
      ),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        size = base_size - 1
      ),
      axis.text.y = element_text(size = base_size - 1),
      axis.title = element_text(size = base_size),
      legend.position = "none"
    )
  
  ## Optinal ylimit
  if (!is.null(y_limits)) {

  limit_df <- do.call(
    rbind,
    lapply(names(y_limits), function(metric) {
      data.frame(
        Metric = metric,
        x = "",
        ymin = y_limits[[metric]][1],
        ymax = y_limits[[metric]][2]
      )
    })
  )

  limit_df$Metric <- factor(
    limit_df$Metric,
    levels = levels(plot_df$Metric)
  )

  p <- p +
    geom_blank(
      data = limit_df,
      aes(x = x, y = ymin),
      inherit.aes = FALSE
    ) +
    geom_blank(
      data = limit_df,
      aes(x = x, y = ymax),
      inherit.aes = FALSE
    )
  }

  ## Optional points
  if (point_size > 0) {
    p <- p +
      geom_jitter(
        width = 0.15,
        size = point_size,
        alpha = 0.3
      )
  }

  return(p)
}