# ==================================================

#' Plot a correlation heatmap from a riboclass object.
#'  
#' Shows the correlation **distance** between samples.
#' @md
#' @inheritParams plot_heatmap
#'
#' @return ComplexHeatmap object
#' @export
#'
#' @examples 
#' data("ribo_toy")
#' plot_heatmap_corr(ribo_toy,"count","condition")
#' 
plot_heatmap_corr <- function(ribo, values_col, color_col=NULL) {
  matrix <- extract_data(ribo, values_col, position_to_rownames = T)
  .plot_heatmap_corr(matrix, ribo[["metadata"]], color_col = color_col)
}
#' Internal function of plot_heatmap_corr.
#'
#' @param cscore_matrix  Sites x Samples C-score matrix (output of extract_data()).
#' @param metadata Metadata of samples in matrix
#' @param color_col Vector of the metadata columns’ name used for coloring samples.
#'
#' @return ComplexHeatmap heatmap
#' 
.plot_heatmap_corr <- function(cscore_matrix, metadata,
                               color_col) {
  
  if(!is.null(color_col)) {
    col <- generate_palette(metadata,color_col)
    column_ha <- ComplexHeatmap::HeatmapAnnotation(df = metadata[color_col], col = col)
  } else {
    column_ha <- NULL
  }
  
  white_red <- colorRamp2::colorRamp2(c(0,1),c("white", "red"))
  
   corr_matrix <- 1 - cor(cscore_matrix,use = "complete.obs")
   dist_cor <- as.dist(corr_matrix)
   
  ComplexHeatmap::Heatmap(corr_matrix,col = white_red,name = "Correlation-based distance",
                          row_title = "Sample",column_title = "Sample", 
                          column_title_side = "bottom",
                          cluster_rows = FALSE, cluster_columns = TRUE,
                          clustering_distance_columns = "manhattan", 
                          clustering_method_columns = "ward.D2",
                          column_split = 3,
                          top_annotation = column_ha,
                          row_names_gp = grid::gpar(fontsize = 6))
  
}
