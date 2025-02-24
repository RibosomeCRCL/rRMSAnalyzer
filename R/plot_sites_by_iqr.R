#' Plot IQR or variance per site / Plot boxplot per site ID. The sites that show an IWR/variance higher than median + 2*MAD are colored in red and are considered as the most variant sites.
#'
#' @param ribo A ribo class object containing the data
#' @param plot Type of plot: "IQR" (default) or "boxplot"
#' @param variance Metric to visualize: "IQR" (default) or "var" (variance)
#' @param ... Additional theme arguments for ggplot
#'
#' @return A ggplot object
#' @export
#'
#' @examples 
#' data("ribo_toy")
#' data("human_methylated")
#' ribo_toy <- rename_rna(ribo_toy)
#' ribo_toy <- annotate_site(ribo_toy,human_methylated)
#' plot_sites_by_IQR(ribo = ribo_toy, plot = "IQR")

plot_sites_by_IQR <- function(ribo = NULL, plot = "IQR", variance = "IQR", ...) {
  
  # Extract C-score matrix from ribo object
  Cscore_matrix <- extract_data(ribo, col = "cscore", position_to_rownames = TRUE, only_annotated = TRUE)
  
  # Compute IQR or variance and order the sites accordingly
  IQR_order_df <- get_IQR(Cscore_matrix, order = variance)
  IQR_order_df <- IQR_order_df[order(IQR_order_df$iqr, decreasing = FALSE), ]
  IQR_order_df$site <- factor(x = IQR_order_df$site, levels = unique(IQR_order_df$site))
  
  # Compute threshold for highlighting significant points (median + 2 * MAD)
  median_2mad <- median(unique(IQR_order_df[,1:2])[,1]) + 2 * mad(unique(IQR_order_df[,1:2])[,1])
  
  # Check if the requested plot type is "IQR" or "var"
  if (tolower(plot) %in% c("iqr", "var")) {
    
    # Scatter plot of variance/IQR per site with color highlighting
    p1 <- ggplot(IQR_order_df, aes_string(x = "site", y = tolower(variance))) +
      geom_point(aes(colour = ifelse(IQR_order_df[,1:2][,1] > median_2mad, "red", "black"))) +
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 45, size = 9, hjust = 1), 
            legend.position = "top", 
            axis.title = element_text(size = 14), 
            axis.text.y = element_text(size = 12)
            ) +
      labs(x = "rRNA 2'Ome sites",
           y = ifelse(tolower(variance) == "iqr", "Interquartile Range (IQR)", "Variance"), size = 20) +
      theme(...) +
      scale_color_identity()
    
    # Density plot of the variance/IQR
    dp <- ggplot(unique(IQR_order_df[,1:2]), aes_string(x = tolower(variance))) + geom_density()
    d <- ggplot_build(dp)$data[[1]]
    
    # Create a density plot alongside the main plot for visualization
    y_density <- cowplot::axis_canvas(p1, axis = "y", coord_flip = TRUE) +
      geom_density(data = unique(IQR_order_df[,1:2]), aes_string(x = tolower(variance))) +
      geom_area(data = subset(d, x >= median_2mad), aes(x = x, y = y), fill = "red", alpha = 0.2) +
      coord_flip()
    
    # Combine the scatter plot and density plot
    combined_plot <- cowplot::insert_yaxis_grob(p1, y_density, position = "right")
    return(plot(combined_plot))
  }
  
  # If "boxplot" is selected, generate a boxplot of C-score per site
  if (tolower(plot) == "boxplot") {
    p1 <- ggplot(IQR_order_df, aes(x = sites.id, y = Cscore)) +
      geom_boxplot() +
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 45, size = 9, hjust = 1), 
            legend.position = "top", 
            axis.title = element_text(size = 14), 
            axis.text.y = element_text(size = 12)
            ) +
      labs(x = "rRNA 2'Ome sites",
           y = "Cscore", size = 20) +
      theme()
    return(p1)
  }
}
