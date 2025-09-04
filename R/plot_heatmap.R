#' Plot heatmap for a riboclass object
#' @description This easy function will let you display an heatmap for any given column (count or c-score). You can add an additionnal layer of information with metadata columns.
#' @param ribo A RiboClass object.
#' @param color_col Vector of the metadata columns’ name used for coloring samples.
#' @param only_annotated Use only annotated sites (default = TRUE).
#' @param title Title to display on the plot. "default" for default title.
#' @param cutree_rows number of clusters the rows are divided into, based on the hierarchical clustering (using cutree).
#' @param cutree_cols number of clusters the columns are divided into, based on the hierarchical clustering (using cutree).
#' @param ... Pheatmap’s parameters
#' @return A ggplot object of a heatmap. The distance used is manhattan and the clustering method is Ward.D2. See ComplexHeatmap doc for more details
#' @export
#'
#' @examples 
#' data("ribo_toy")
#' data("human_methylated") 
#' ribo_toy <- rename_rna(ribo_toy)
#' ribo_toy <- annotate_site(ribo_toy,
#'                                 annot = human_methylated,
#'                                 anno_rna = "rRNA",
#'                                 anno_pos = "Position",
#'                                 anno_value = "Nomenclature")
#' plot_heatmap(ribo_toy,  color_col = c("run","condition"), only_annotated=TRUE)
plot_heatmap <- function(ribo, color_col = NULL, only_annotated=FALSE, title,
                         cutree_rows=4, cutree_cols=2, ...) {
  
  check_metadata(ribo,color_col)
  matrix <- extract_data(ribo, "cscore", position_to_rownames = TRUE,
                         only_annotated = only_annotated)
  if (nrow(matrix) == 0) message("Empty matrix : Are your data annotated ?")
  
  .plot_heatmap(matrix, ribo[["metadata"]], color_col = color_col,
                most_variant = FALSE, title = title, cutree_rows = cutree_rows,
                cutree_cols = cutree_cols,...)
}


#' Internal function to plot heatmap
#'
#' @param cscore_matrix Sites x Samples C-score matrix (output of extract_data()).
#' @param metadata metadata for samples in cscore_matrix
#' @inheritParams plot_heatmap
#' @param most_variant select only the most variant positions (cannot be used from plot_heatmap())
#'
#' @return ComplexHeatmap heatmap
#' @keywords internal
.plot_heatmap <- function(cscore_matrix = NULL, metadata = NULL,
                                  color_col = NULL, most_variant = FALSE,
                                  title="default", cutree_rows,
                                  cutree_cols, ...) {

  heat_colors <- grDevices::hcl.colors(7,"inferno")
  
  
  if(!is.null(color_col)) {
    col <- generate_palette(metadata,color_col)
    column_ha <- ComplexHeatmap::HeatmapAnnotation(df = metadata[color_col], col = col, na_col = "red")
  } else {
    column_ha <- NULL
  }
  
  cscore_matrix <- stats::na.omit(cscore_matrix)
  cscore_matrix <- as.matrix(cscore_matrix)
ComplexHeatmap::Heatmap(cscore_matrix,col = heat_colors,name = "C-score",
                          row_title = "rRNA 2'Ome sites",column_title = "Samples", 
                          column_title_side = "bottom",
                          cluster_rows = TRUE, cluster_columns = TRUE,
                          clustering_distance_columns = "manhattan", 
                          clustering_distance_rows = "manhattan",
                          clustering_method_columns = "ward.D2",
                          clustering_method_rows = "ward.D2",
                          row_split = cutree_rows, column_split = cutree_cols,
                          top_annotation = column_ha,
                          row_names_gp = grid::gpar(fontsize = 6) 
                          )
}