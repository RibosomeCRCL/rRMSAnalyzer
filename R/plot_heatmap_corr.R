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
plot_heatmap_corr <- function(ribo, values_col, color_col) {
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
.plot_heatmap_corr <- function(cscore_matrix = NULL, metadata = NULL,
                               color_col = NULL) {
  
  # rownames of metadata
  rownames(metadata) <- metadata[, "samplename"]
  
  # all character columns to factor columns
  metadata[sapply(metadata, is.character)] <- lapply(
    metadata[sapply(metadata, is.character)],
    as.factor
  )

  corr_matrix <- 1 - cor(cscore_matrix,use = "complete.obs")
  dist_cor <- as.dist(corr_matrix)
  metadata_1 <- NA
  if(!is.null(color_col)) {
    metadata_1 <- metadata[color_col]
    metadata_1 <- data.frame(metadata_1)
    
  }
  else {
    metadata_1 <- NULL
  }
  
  white_red <- colorRampPalette(c("white", "red"), interpolate = "linear")(100)

    htmap <- ComplexHeatmap::pheatmap(corr_matrix,
                                clustering_method = "ward.D2",
                                cluster_rows = FALSE,
                                clustering_distance_cols = dist_cor,
                                clustering_distance_rows = dist_cor,
                                color = white_red,
                                breaks = seq(0, 1, by = 0.01),
                                annotation_col = metadata_1,
                                main = "Correlation-based distance heatmap"
    )
  
  return(htmap)
}
