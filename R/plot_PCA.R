#' Title
#'
#' @param ribo 
#' @param col_to_plot 
#' @param order_by_col 
#' @param col_for_color 
#'
#' @return
#' @export
#'
#' @examples
plot_PCA <- function(ribo,col_to_plot,order_by_col,col_for_color = NULL) {
  
  pca_matrix <- aggregate_samples_by_col(ribo[["raw_counts"]],col_to_plot,position_to_rownames = T)
  
  pca_calculated <- .calculate_pca(pca_matrix, ribo[["metadata"]],order_by_col)
  
  return(.plot_pca(pca_calculated,ribo[["metadata"]],col_for_color))
  
}


#' Title
#'
#' @param cscore.matrix 
#' @param metadata 
#' @param order.by.col 
#'
#' @return
#' @export
#'
#' @examples
.calculate_pca <- function(cscore.matrix = NULL, metadata = NULL, order.by.col = NULL) {
  pca.res <- dudi.pca(t(cscore.matrix[complete.cases(cscore.matrix), match(metadata[,order.by.col], colnames(cscore.matrix))]), 
                      scannf = F, 
                      nf = 5)
  
  return(pca.res)
}


#' Title
#'
#' @param dudi.pca 
#' @param metadata 
#' @param col.by.col 
#'
#' @return
#' @export
#'
#' @examples
.plot_pca <- function(dudi.pca = NULL, metadata = NULL, col.by.col = NULL) {
  plot.pca <- factoextra::fviz_pca_ind(dudi.pca, 
                           title = paste("PCA of Cscore for", as.character(ncol(dudi.pca$tab)), "sites."),
                           repel = T, 
                           habillage = metadata[,col.by.col],
                           pointsize = 2, 
                           labelsize = 4) + theme(text = element_text(size = 18)) + labs(subtitle = paste(length(metadata[,"filename"]),"samples"))
  
  return(plot.pca)
}
