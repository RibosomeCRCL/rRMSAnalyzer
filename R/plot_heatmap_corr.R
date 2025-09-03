#' Plot a correlation heatmap from a riboclass object
#'  
#' Shows the correlation distance between samples.
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
#' Internal function of plot_heatmap_corr
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
  
  if(!is.null(color_col)) { # if color is given
    col <- generate_palette(metadata,color_col) # generate a color palette for the heatmap
    column_ha <- ComplexHeatmap::HeatmapAnnotation(df = metadata[color_col], col = col) # annotation creation
  } else {
    column_ha <- NULL
  }
    corr_matrix <- stats::cor(cscore_matrix,use = "complete.obs") # calculate correlation between lines or columns matrix, correlation is calculated ignoring missing values (NA)
    pearson_color <- colorRamp2::colorRamp2(c(0,0.5,1),c("red","white", "blue")) # 0 = red, 0.5 = white, 1 = blue
   
 ht <- ComplexHeatmap::Heatmap(corr_matrix,col = pearson_color,name = "Pearson correlation",
                          row_title = "Sample",column_title = "Pearson correlation between samples", 
                          column_title_side = "top", # column title position on top
                          cluster_rows = FALSE, cluster_columns = TRUE, # no clustering of lines (samples stay in original order), column clustering (correlated samples will be grouped
                          clustering_distance_columns = "pearson", # use of pearson distance to mesure similarity between samples
                          heatmap_legend_param = list(
                            legend_direction = "horizontal" ),
                          clustering_method_columns = "ward.D2", # use of clustering method Ward.D2 which group samples minimizing intra-group variance
                          column_split = 3, # divide colums in 3 groups after clustering
                          top_annotation = column_ha, # add annotations defined in the loop on the top
                          row_names_gp = grid::gpar(fontsize = 6))
  
 ComplexHeatmap::draw(ht, heatmap_legend_side = "bottom") # display legend in a horizontal manner
}
