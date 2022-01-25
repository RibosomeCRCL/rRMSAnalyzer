#' Title
#'
#' @param ribo 
#' @param col_to_plot 
#' @param order_by_col 
#' @param col_for_color
#' @param axis
#'
#' @return
#' @export
#'
#' @examples
plot_PCA <- function(ribo,col_to_plot,order_by_col,col_for_color = NULL, axis = c(1,2)) {
  
  # Faut faire des controles de chaque params
  
  pca_matrix <- aggregate_samples_by_col(ribo[["counts"]],col_to_plot,position_to_rownames = T)
  
  pca_calculated <- .calculate_pca(pca_matrix, ribo[["metadata"]],order_by_col)
  
  return(.plot_pca(pca_calculated,ribo[["metadata"]],col_for_color, axis = axis))
  
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
.plot_pca <- function(dudi.pca = NULL, metadata = NULL, col.by.col = NULL, axis = axis) {
  plot.pca <- factoextra::fviz_pca_ind(dudi.pca, 
                           title = paste("PCA of Cscore for", as.character(ncol(dudi.pca$tab)), "sites."),
                           repel = TRUE, 
                           habillage = metadata[,col.by.col],
                           pointsize = 2, 
                           labelsize = 4,
                           axes = axis) + 
    theme_bw()+
    theme(text = element_text(size = 16)) + 
    labs(subtitle = paste(as.character(nrow(dudi.pca$tab)),"samples")) + 
    ylab(paste("PC", axis[2],":", round(dudi.pca$eig[axis[2]]/sum(dudi.pca$eig) * 100, 1), "%")) +
    xlab(paste("PC", axis[1],":", round(dudi.pca$eig[axis[1]]/sum(dudi.pca$eig) * 100, 1), "%"))
    
  
  return(plot.pca)
}
