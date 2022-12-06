#' Principal component analysis of a RiboClass object
#'
#' 
#' @param ribo A RiboClass object.
#' @param color_col Name of the column in the metadata used for coloring PCA inds.
#' @param axes Two-element vector indicating which pair of principal components you want to show.
#' @param only_annotated If TRUE, use only annotated sites to plot PCA.
#' @param title tile to display on the plot. "default" for default title.
#' @param subtitle subtitle to display on the plot. "samples" for number of samples. "none" for no subtitle.
#' @param pca_object_only Return directly the full dudi.pca object, without generating the plot.
#' @return A ggplot or a dudi.pca object if pca_object_only is set to True.
#' @export
#'
#' @examples 
#' data("ribo_toy")
#' plot_PCA(ribo_toy,"run")
plot_PCA <- function(ribo, color_col = NULL, axes = c(1,2), only_annotated = FALSE, title = "default", subtitle = "samples", pca_object_only = FALSE) {
  
  if (is.null(ribo)) {stop("MISSING parameter: please provide a RiboClass!")}
  if (!inherits(ribo, "RiboClass")) {stop("ribo argument is not a RiboClass!")}
  
  
  if (isFALSE(ribo$has_cscore)) {stop("You should calculate Cscores first using calculate_score funciton")}

  
  pca_matrix <- extract_data(ribo,"cscore",position_to_rownames = TRUE)
  
  if(only_annotated) pca_matrix <- pca_matrix[which(!is.na(ribo[["data"]][[1]][["site"]])),]
    
  
  pca_calculated <- .calculate_pca(pca_matrix)
  
  if(pca_object_only) return(pca_calculated)
  
  facto_pca <- .plot_pca(pca_calculated,ribo[["metadata"]],color_col, axes = axes, title = title, subtitle = subtitle)
  
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
.plot_pca <- function(dudi.pca = NULL, metadata = NULL, color_col = NULL, axes = axes, title="default", subtitle = "samples") {
  # col.by.com argument can be NULL if no metadata is given by the user
  if(is.null(color_col)) {
    color_column <- "none"
  }
  else {
    color_column <- metadata[,color_col]
  }
  # Plot title
  if(title == "default") {
    plot_title <- paste("PCA of Cscore for", as.character(ncol(dudi.pca$tab)), "sites")
  }
  else {
    plot_title <- title
  }
  # Plot subtitle
  
  if(subtitle == "samples") {
   plot_subtitle <- paste(as.character(nrow(dudi.pca$tab)),"samples")
  }
  else if(subtitle == "none") {
    plot_subtitle <- ggplot2::waiver()
  }
  else {
    plot_subtitle <- subtitle
  }
  
  plot.pca <- factoextra::fviz_pca_ind(dudi.pca, 
                                       title = plot_title ,
                                       repel = TRUE, 
                                       habillage = color_column,
                                       pointsize = 2, 
                                       labelsize = 4,
                                       axes = axes) + 
    ggplot2::theme_bw()+
    ggplot2::theme(text = ggplot2::element_text(size = 16)) + 
    ggplot2::labs(subtitle = plot_subtitle) + 
    ggplot2::ylab(paste("PC", axes[2],":", round(dudi.pca$eig[axes[2]]/sum(dudi.pca$eig) * 100, 1), "%")) +
    ggplot2::xlab(paste("PC", axes[1],":", round(dudi.pca$eig[axes[1]]/sum(dudi.pca$eig) * 100, 1), "%"))
  # TODO PCA legend should be according to metadata
  
  
  return(plot.pca)
}
