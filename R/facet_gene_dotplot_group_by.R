facet_gene_dotplot_group_by <- function(
    cds,
    markers,
    group_cells_by = "clusters",      # y-axis grouping
    facet_by = "Timepoint",           # x-axis grouping
    reduction_method = "UMAP",
    norm_method = c("log", "size_only"),
    lower_threshold = 0,
    max.size = 10,
    axis_order = c('group_marker', 'marker_group'),
    flip_percentage_mean = FALSE,
    pseudocount = 1,
    scale_max = 3,
    scale_min = -3
) {
  
  # Get gene IDs for markers
  gene_ids <- fData(cds) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "rowname") %>%
    dplyr::filter(rowname %in% markers | gene_short_name %in% markers) %>%
    dplyr::pull(rowname)
  
  # Expression matrix
  exprs_mat <- t(as.matrix(exprs(cds)[gene_ids, ]))
  exprs_mat <- reshape2::melt(exprs_mat)
  colnames(exprs_mat) <- c('Cell', 'Gene', 'Expression')
  exprs_mat$Gene <- as.character(exprs_mat$Gene)
  
  # Extract grouping columns from colData
  cell_group <- colData(cds)[, group_cells_by]
  names(cell_group) <- colnames(cds)
  exprs_mat$Group <- cell_group[exprs_mat$Cell]
  
  facet_group <- colData(cds)[, facet_by]
  names(facet_group) <- colnames(cds)
  exprs_mat$Facet <- facet_group[exprs_mat$Cell]
  
  # Remove NAs
  exprs_mat <- exprs_mat %>% dplyr::filter(!is.na(Group) & !is.na(Facet))
  
  # Compute summary: mean log expression + fraction expressing
  ExpVal <- exprs_mat %>%
    dplyr::group_by(Group, Gene, Facet) %>%
    dplyr::summarize(
      mean = mean(log(Expression + pseudocount)),
      percentage = sum(Expression > lower_threshold) / length(Expression),
      .groups = "drop"
    )
  
  # Clamp mean values
  ExpVal$mean <- pmin(pmax(ExpVal$mean, scale_min), scale_max)
  
  # Convert Gene IDs to gene_short_name
  ExpVal$Gene <- fData(cds)[ExpVal$Gene, 'gene_short_name']
  
  # Plot
  g <- ggplot(ExpVal, aes(y = Group, x = Facet)) +
    geom_point(aes(colour = mean, size = percentage)) +
    viridis::scale_color_viridis(name = paste0('log(mean + ', pseudocount, ')')) +
    scale_size(name = 'percentage', range = c(0, max.size)) +
    facet_wrap(~Gene, nrow =1) +
    theme_bw() +
    xlab(facet_by) + ylab(group_cells_by) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(g)
}