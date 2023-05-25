#' Principal component analysis of a RiboClass object
#'
#' 
#' @param ribo A RiboClass object.
#' @param color_col Name of the column in the metadata used for coloring samples.
#' @param axes Two-element vector indicating which pair of principal components 
#' to show.
#' @param only_annotated If TRUE, use only annotated sites to plot PCA.
#' @param title Title to display on the plot. 'default' for default title.
#' @param subtitle Subtitle to display on the plot. 'samples' for number of
#' samples. 'none' for no subtitle.
#' @param object_only Return directly the full dudi.pca object, without
#' generating the plot.
#' @return A ggplot or a dudi.pca object if object_only is set to True.
#' @export
#'
#' @examples 
#' data('ribo_toy')
#' plot_pca(ribo_toy,'run')
plot_pca <- function(ribo, color_col = NULL, axes = c(1, 2),
                     only_annotated = FALSE, title = "default",
                     subtitle = "samples", object_only = FALSE) {
  
  if (missing(ribo)) {
    cli::cli_abort(c(
      "Missing argument : {.var ribo}",
      "i" = "{.var ribo} only accepts a RiboClass object."
    ))
  }
  if (!inherits(ribo, "RiboClass")) {
    cli::cli_abort(c(
      "{.var ribo} must be a RiboClass object",
      "x" = "You've supplied a {.cls {class(ribo)}}."
    ))
  }
  if (isFALSE(ribo$has_cscore)) {
    cli::cli_abort(c(
      "No C-score found in the RiboClass supplied in {.var ribo}!",
      "i" = "You can compute C-scores using compute_cscore function."
      ))
  }
  
  if (!is.null(color_col))
    check_metadata(ribo, color_col)
  
  pca_matrix <- extract_data(ribo, "cscore", position_to_rownames = TRUE,
                             only_annotated = only_annotated)
  
  pca_calculated <- .calculate_pca(pca_matrix)
  
  if (object_only)
    return(pca_calculated)
  
  facto_pca <- .plot_pca(pca_calculated, ribo[["metadata"]], color_col,
                         axes = axes, title = title, subtitle = subtitle)
  return(facto_pca)
}

#' Compute PCA from a c-score matrix
#'
#' @keywords internal
#'
#' @param cscore.matrix matrix of c-score extracted from a RiboClass with ' \code{\link{extract_data}}
#' @return dudi.pca object
#'
.calculate_pca <- function(cscore.matrix = NULL) {
  pca.res <- ade4::dudi.pca(t(cscore.matrix[stats::complete.cases(cscore.matrix),]), 
                            scannf = FALSE, 
                            nf = 5)
  
  return(pca.res)
}

#' Plot a dudi.pca object using factoextra
#' 
#' @keywords internal
#'
#' @param dudi.pca a dudi.pca object generated with Ade4
#' @param metadata metadata table from RiboClass
#' @inheritParams plot_pca
#'
#' @return a ggplot
#'
.plot_pca <- function(dudi.pca = NULL, metadata = NULL,
                      color_col = NULL, axes = axes, title = "default",
                      subtitle = "samples") {

  if (is.null(color_col)) {
    color_column <- "none"
  } else {
    color_column <- metadata[, color_col]
  }
  # Plot title
  if (title == "default") {
    plot_title <- "Principal Component Analysis from C-score data"
  } else {
    plot_title <- title
  }
  # Plot subtitle
  
  if (subtitle == "samples") {
    plot_subtitle <- paste(as.character(nrow(dudi.pca$tab)),
                           "samples and", as.character(ncol(dudi.pca$tab)),
                           "positions")
  } else if (subtitle == "none") {
    plot_subtitle <- ggplot2::waiver()
  } else {
    plot_subtitle <- subtitle
  }
  
  plot.pca <- factoextra::fviz_pca_ind(dudi.pca,
                                       title = plot_title, repel = TRUE,
                                       habillage = color_column, pointsize = 2,
                                       labelsize = 4, axes = axes) +
    ggplot2::theme_bw() +
    ggplot2::theme(text = ggplot2::element_text(size = 12)) +
    ggplot2::labs(subtitle = plot_subtitle) +
    ggplot2::ylab(paste("PC", axes[2], ":",
                        round(dudi.pca$eig[axes[2]]/sum(dudi.pca$eig) * 100, 1),
                        "%")) +
    ggplot2::xlab(paste("PC", axes[1], ":",
                        round(dudi.pca$eig[axes[1]]/sum(dudi.pca$eig) * 100, 1),
                        "%"))
  return(plot.pca)
}
