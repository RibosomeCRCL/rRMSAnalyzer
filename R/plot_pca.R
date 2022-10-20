#' Principal component analysis of a RiboClass object
#'
#' @param ribo a RiboClass object
#' @param color_col column in the metadata used for coloring the PCA.
#' @param axes two-element vector indicating which pair of principal components you want to show.
#' @param only_annotated Use only annotated sites to plot PCA.
#' @param pca_object_only Return directly the full dudi.pca object, without generating the plot.
#' @return a ggplot or a dudi.pca object if pca_object_only is set to True.
#' @export
#'
#' @examples 
#' data("ribo_toy")
#' plot_PCA(ribo_toy,"run")
plot_PCA <- function(ribo, color_col = NULL, axes = c(1,2), only_annotated = FALSE, pca_object_only = FALSE) {
  
  if (is.null(ribo)) {stop("MISSING parameter: please provide a RiboClass!")}
  if (is(ribo) != "RiboClass") {stop("ribo argument is not a RiboClass!")}
  
  
  if (isFALSE(ribo$has_cscore)) {stop("You should calculate Cscores first using calculate_score funciton")}

  
  pca_matrix <- extract_data(ribo,"cscore",position_to_rownames = TRUE)
  
  if(only_annotated) pca_matrix <- pca_matrix[which(!is.na(ribo[["data"]][[1]][["site"]])),]
    
  
  pca_calculated <- .calculate_pca(pca_matrix)
  
  if(pca_object_only) return(pca_calculated)
  
  facto_pca <- .plot_pca(pca_calculated,ribo[["metadata"]],color_col, axes = axes)
  
  return(facto_pca)
  
}




#' Plot a dudi.pca object using factoextra
#' 
#' @keywords internal
#'
#' @param dudi.pca a dudi.pca object generated with Ade4
#' @param metadata metadata table from RiboClass
#' @inheritParams plot_PCA
#'
#' @return a ggplot
#'
.plot_pca <- function(dudi.pca = NULL, metadata = NULL, color_col = NULL, axes = axes) {
  # col.by.com argument can be NULL if no metadata is given by the user
  if(is.null(color_col)) {
    color_col <- "none"
  }
  else {
    color_col <- metadata[,color_col]
  }
  plot.pca <- factoextra::fviz_pca_ind(dudi.pca, 
                                       title = paste("PCA of Cscore for", as.character(ncol(dudi.pca$tab)), "sites"),
                                       repel = TRUE, 
                                       habillage = color_col,
                                       pointsize = 2, 
                                       labelsize = 4,
                                       axes = axes) + 
    ggplot2::theme_bw()+
    ggplot2::theme(text = ggplot2::element_text(size = 16)) + 
    ggplot2::labs(subtitle = paste(as.character(nrow(dudi.pca$tab)),"samples")) + 
    ggplot2::ylab(paste("PC", axes[2],":", round(dudi.pca$eig[axes[2]]/sum(dudi.pca$eig) * 100, 1), "%")) +
    ggplot2::xlab(paste("PC", axes[1],":", round(dudi.pca$eig[axes[1]]/sum(dudi.pca$eig) * 100, 1), "%"))
  # TODO PCA legend should be according to metadata
  
  
  return(plot.pca)
}
