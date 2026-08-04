plot_go_dotplot <- function(df, total_genes = NULL, top_n = 20) {
  
  # Check required columns
  required_cols <- c("p_value", "term_id", "term_name", "intersection_size")
  if(!all(required_cols %in% colnames(df))) {
    stop("Input dataframe must contain columns: p_value, term_id, term_name, intersection_size")
  }
  
  # Compute gene ratio
  if(is.null(total_genes)) {
    stop("Please provide total_genes (number of DEGs or total background genes)")
  }
  df <- df %>%
    mutate(gene_ratio = intersection_size / total_genes)
  
  # Optional: take top N terms by significance
  df <- df %>%
    arrange(p_value) %>%
    slice_head(n = top_n)
  
  # Reorder terms for plotting
  df$term_name <- factor(df$term_name, levels = rev(df$term_name))
  
  # Create dotplot
  p <- ggplot(df, aes(x = gene_ratio, y = term_name)) +
  geom_point(aes(size = intersection_size, color = -log10(p_value))) +
  scale_color_gradient(low = "lightblue", high = "darkblue", name = "-log10(p_value)") +
  scale_size_continuous(name = "Intersection size") +
  labs(x = "Gene ratio", y = "GO term", title = "GO Term Enrichment Dotplot") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )
  
  return(p)
}