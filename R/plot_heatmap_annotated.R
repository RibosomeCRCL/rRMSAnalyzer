#' Plot a correlation heatmap from a RiboClass object.
#' @description
#' Shows the correlation distance between samples.
#'
#' @param ribo A RiboClass object.
#' @param qcdata a dataframe with metadata and outliers information
#'
#' @returns ComplexHeatmap object
#' @export
#'
#' @examples 
#' data("ribo_toy")
#' qcdata <- ribo_toy[[2]]
#' ribo_matrix <- extract_data(ribo_toy, "count", position_to_rownames = TRUE)
#' medcov <- apply(ribo_matrix,2,function(x) median(x,na.rm=TRUE)) 
#' qcdata <- cbind(qcdata,median_coverage=medcov[rownames(qcdata)]) 
#' qcdata$coverage_quality <- "pass" 
#' qcdata$coverage_quality[qcdata$median_coverage < 100] <- "warning" 
#' plot_rle <- rRMSAnalyzer::plot_rle(ribo_toy,show_outlier = TRUE)
#' mad <- plot_rle$plot_env$mad  
#' rle_median <- as.vector(plot_rle$plot_env$rle_grouped$median)
#' keys <- as.vector(plot_rle$plot_env$rle_grouped$key)
#' outlier_table <- data.frame(key = keys, rle_median = rle_median)
#' outlier <- data.frame()
#' for (i in 1:nrow(outlier_table)) {
#'   if(outlier_table$rle_median[i] < mad) {
#'     outlier <- rbind(outlier, outlier_table[i,])}}
#' qcdata$rle_median <- outlier_table$rle_median[match(qcdata$samplename, outlier_table$key)]
#' qcdata$rle_median_quality <- "pass" 
#' qcdata$rle_median_quality[qcdata$rle_median < mad] <- "warning" 
#' qcdata$total_outliers <- rowSums(qcdata[, c("coverage_quality", "rle_median_quality")] == "warning")
#' ribo_toy$metadata$outlier_level <- qcdata$total_outliers
#' qcdata$total_outliers <- as.character(qcdata$total_outliers)
#' plot_heatmap_annotated(ribo_toy, qcdata)

plot_heatmap_annotated <- function(ribo = ribo, qcdata = qcdata) {
  
  corr_matrix <- extract_data(ribo, "count")
  corr_matrix <- corr_matrix[, -1]
  corr_matrix <- cor(corr_matrix, use = "complete.obs")
  pearson_color <- colorRamp2::colorRamp2(c(0,0.5,1),c("red","white", "blue"))
  
  # Préparer les annotations
  row_annotation <- ComplexHeatmap::HeatmapAnnotation(
    outlier_coverage_distribution = qcdata$coverage_quality,
    outlier_relative_log_coverage = qcdata$rle_median_quality,
    col = list(
      outlier_coverage_distribution = c("pass" = "darkgrey", "warning" = "red"),
      outlier_relative_log_coverage = c("pass" = "darkgrey", "warning" = "red")
    ),
    annotation_legend_param = list(
      title_gp = grid::gpar(fontsize = 16), # title size
      labels_gp = grid::gpar(fontsize = 14), # tag title
      grid_width = unit(8, "mm"),     # square width
      grid_height = unit(8, "mm")     # square height
    )
  )
  
  # create a vector for sample color by level of outlier
  column_colors <- ifelse(qcdata$total_outliers > 0, "red", "black")
  
  # identifying warning sample
  column_colors <- ifelse(
    qcdata$coverage_quality == "warning" | qcdata$rle_median_quality == "warning",
    "red",  # "Warning" in red
    "black" # other in black
  )
  
  # heatmap plot with name of coloured columns
  x <- ComplexHeatmap::Heatmap(
    corr_matrix, col = pearson_color,
    top_annotation = row_annotation,
    name = "Pearson correlation",
    column_title = "Pearson correlation between samples", 
    column_title_side = "top",
    cluster_rows = FALSE, cluster_columns = TRUE,
    clustering_distance_columns = "pearson", 
    heatmap_legend_param = list(
      legend_direction = "horizontal", # horizontal legend
      title_gp = grid::gpar(fontsize = 16),  # legend title size
      labels_gp = grid::gpar(fontsize = 14), # legend tags title
      grid_height = unit(10, "mm"),    # legend square height 
      grid_width = unit(10, "mm")      # legend square width
    ),
    clustering_method_columns = "ward.D2",
    column_split = 3,
    row_names_gp = grid::gpar(fontsize = 16),
    rect_gp = grid::gpar(col = "black", lwd = 0.2),
    column_dend_side = "bottom",
    column_dend_height = unit(1, "cm"),
    column_dend_gp = grid::gpar(col = "red"),
    column_names_gp = grid::gpar(col = column_colors, fontsize = 16) # column name color
  )
  
  ComplexHeatmap::draw(x,
       heatmap_legend_side = "bottom",  # legend position
       annotation_legend_side = "bottom" # for legend annotation
  ) 
}