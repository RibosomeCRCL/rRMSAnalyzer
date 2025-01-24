# ==================================================

#' Plot a correlation heatmap from a riboclass object.
#'  
#' Shows the correlation **distance** between samples.
#' @md
#' @inheritParams plot_heatmap
#' @param values_col Name of the column containing the value (either count or cscore).
#' @return ComplexHeatmap object
#' @export
#'
#' @examples 
#' data("ribo_toy")
#' plot_heatmap_corr(ribo_toy,"count","condition")
#' 
plot_heatmap_corr <- function(ribo, values_col, color_col=NULL) {
  matrix <- extract_data(ribo, values_col, position_to_rownames = TRUE)
  if(!is.null(color_col)) check_metadata(ribo,color_col)
  .plot_heatmap_corr(matrix, ribo[["metadata"]], color_col = color_col)
}
#' Internal function of plot_heatmap_corr.
#'
#' @param cscore_matrix  Sites x Samples C-score matrix (output of extract_data()).
#' @param metadata Metadata of samples in matrix
#' @param color_col Vector of the metadata columns’ name used for coloring samples.
#'
#' @return ComplexHeatmap heatmap
#' @keywords internal
#' 
.plot_heatmap_corr <- function(cscore_matrix, metadata,
                               color_col) {
  
  if(!is.null(color_col)) {
    col <- generate_palette(metadata,color_col)
    column_ha <- ComplexHeatmap::HeatmapAnnotation(df = metadata[color_col], col = col)
  } else {
    column_ha <- NULL
  }
    corr_matrix <- stats::cor(cscore_matrix,use = "complete.obs")
    pearson_color <- colorRamp2::colorRamp2(c(0,0.5,1),c("red","white", "blue"))
   
 ht <- ComplexHeatmap::Heatmap(corr_matrix,col = pearson_color,name = "Pearson correlation",
                          row_title = "Sample",column_title = "Pearson correlation between samples", 
                          column_title_side = "top",
                          cluster_rows = FALSE, cluster_columns = TRUE,
                          clustering_distance_columns = "pearson", 
                          heatmap_legend_param = list(
                            legend_direction = "horizontal" ),
                          clustering_method_columns = "ward.D2",
                          column_split = 3,
                          top_annotation = column_ha,
                          row_names_gp = grid::gpar(fontsize = 6))
  
 ComplexHeatmap::draw(ht, heatmap_legend_side = "bottom")
  
}
