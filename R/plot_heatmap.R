#' plot heatmap for a ribo
#'
#' @param ribo 
#' @param col_to_plot 
#' @param order_by_col 
#' @param cols_for_annotation 
#' @param most_variant 
#'
#' @return
#' @export
#'
#' @examples
plot_heatmap <- function(ribo,col_to_plot,cols_for_annotation, order_by_col,most_variant=F) {
  
  matrix <- aggregate_samples_by_col(ribo[["raw_counts"]],col_to_plot,position_to_rownames = T)
  .plot_heatmap(matrix,ribo[["metadata"]],order_by_col,cols_for_annotation = cols_for_annotation, most_variant = most_variant)
  
}

#' Title
#'
#' @param cscore_matrix 
#' @param annotation_samples 
#' @param order_by_col 
#' @param cols_for_annotation 
#' @param most_variant 
#'
#' @return
#' @export
#'
#' @examples
.plot_heatmap <- function(cscore_matrix = NULL, annotation_samples = NULL,cols_for_annotation = NULL, order_by_col = NULL, most_variant = F) {
  
  annotation_samples <- annotation_samples[,c(order_by_col, cols_for_annotation)]
  # rownames of annotation_samples 
  rownames(annotation_samples) <- annotation_samples[,order_by_col]
  
  
  # all character columns to factor columns
  annotation_samples[sapply(annotation_samples, is.character)] <- lapply(annotation_samples[sapply(annotation_samples, is.character)], 
                                                                         as.factor)
  
  annotation_samples.1 <- annotation_samples[cols_for_annotation]
  # library(RColorBrewer)
  # list.cols <- list()
  #  for (i in seq_along(colnames(annotation_samples.1))) {
  #    
  #  list.cols[[colnames(annotation_samples.1)[i]]] <- setNames(brewer.pal(nlevels(annotation_samples.1[, c(i)]), paste("Set", i, sep = "")), levels(annotation_samples.1[,c(i)]))
  #    
  #  }
  
  if(most_variant) {
    cscore_matrix <- get_most_variant(cscore_matrix)
  }

  
  htmap <- pheatmap::pheatmap(cscore_matrix[complete.cases(cscore_matrix),match(annotation_samples[,order_by_col], colnames(cscore_matrix))] ,
                              clustering_method = "ward.D2" ,
                              clustering_distance_cols = "manhattan", 
                              clustering_distance_rows = "manhattan", 
                              cutree_cols = 2, 
                              cutree_rows = 4, 
                              annotation_col = data.frame(annotation_samples.1),main = "heatmap")
  return(htmap)
  
}

# ==================================================

#' plot a correlation heatmap from a riboclass object
#'
#' @param ribo 
#' @param col_to_plot 
#' @param cols_for_annotation 
#' @param order_by_col 
#'
#' @return
#' @export
#'
#' @examples
plot_heatmap_corr <- function(ribo,col_to_plot,cols_for_annotation, order_by_col,use_triangle=F) {
  
  matrix <- aggregate_samples_by_col(ribo[["raw_counts"]],col_to_plot,position_to_rownames = T)
  .plot_heatmap_corr(matrix,ribo[["metadata"]],order.by.col=order_by_col,cols_for_annotation = cols_for_annotation,use_triangle)
  
}
#' Internal function of plot_heatmap_corr. 
#'
#' @param cscore.matrix 
#' @param annotation.samples 
#' @param cols_for_annotation 
#' @param order.by.col 
#' @param use_triangle 
#'
#' @return
#' @export
#'
#' @examples
.plot_heatmap_corr <- function(cscore.matrix = NULL, annotation.samples = NULL,cols_for_annotation = NULL ,order.by.col = NULL,use_triangle=F) {
  
  # rownames of annotation.samples 
  rownames(annotation.samples) <- annotation.samples[,order.by.col]
  
  
  # all character columns to factor columns
  annotation.samples[sapply(annotation.samples, is.character)] <- lapply(annotation.samples[sapply(annotation.samples, is.character)], 
                                                                         as.factor)
  corr_matrix <- 1- cor(cscore.matrix)
  dist.cor <- as.dist(corr_matrix)
  annotation.samples.1 <- annotation.samples[cols_for_annotation]
  # library(RColorBrewer)
  # list.cols <- list()
  # for (i in seq_along(colnames(annotation.samples.1))) {
  #   
  #   list.cols[[colnames(annotation.samples.1)[i]]] <- setNames(brewer.pal(nlevels(annotation.samples.1[, c(i)]), paste("Set", i, sep = "")), levels(annotation.samples.1[,c(i)]))
  #   
  # }
  white_red<-colorRampPalette(c("white","red"),interpolate = "linear")(100)
  
  if(use_triangle) {
  corr_matrix[lower.tri(corr_matrix)] <- NA
  
  htmap <- pheatmap::pheatmap(corr_matrix[,match(annotation.samples[,order.by.col], colnames(corr_matrix))] ,
                              cluster_cols = F,
                              cluster_rows = F,
                              color = white_red,
                              breaks = seq(0, 1, by = 0.01),
                              annotation_col = data.frame(annotation.samples.1),main = "correlation heatmap")
  }
  else  {
    htmap <- pheatmap::pheatmap(corr_matrix[,match(annotation.samples[,order.by.col], colnames(corr_matrix))] ,
                                clustering_method = "ward.D2" ,
                                clustering_distance_cols = dist.cor, 
                                clustering_distance_rows = dist.cor, 
                                cutree_cols = 2, 
                                cutree_rows = 4, 
                                color = white_red,
                                breaks = seq(0, 1, by = 0.01),
                                annotation_col = data.frame(annotation.samples.1),main = "correlation heatmap")
    
  }
  return(htmap)
  
}

get_most_variant <- function(df = NULL, column_or_row = "row", n = 20,i_want_most_var = TRUE) {
  
  df <- as.data.frame(df)
  
  if(column_or_row == "row") {column_or_row <- 1}
  
  if(column_or_row == "column") {column_or_row <- 2}
  
  if(is.null(df)) {stop("No matix in entry")}
  
  var <- apply(df, column_or_row, var)
  
  df[order(var, decreasing = i_want_most_var)[1:n],]
  
}
