#' Correspondence analysis of a RiboClass' counts.
#'
#' @param ribo A RiboClass object.
#' @param color_col Name of the column in the metadata used for coloring samples.
#' @param axes Two-element vector indicating which pair of COA components
#' to show.
#' @param only_annotated If TRUE, use only annotated sites to plot COA.
#' @param title Title to display on the plot. 'default' for default title.
#' @param subtitle Subtitle to display on the plot. 'samples' for number of
#' samples. 'none' for no subtitle.
#' @param object_only Return directly the full dudi.coa object, without
#' generating the plot.
#' @return A ggplot or a dudi.coa object if object_only is set to True.
#' @export
#'
#' @examples 
#' data('ribo_toy')
#' plot_coa(ribo = ribo_toy, color_col = 'condition')
plot_coa <- function(ribo, color_col = NULL, axes = c(1, 2),
                     only_annotated = FALSE, title = "default",
                     subtitle = "default", object_only = FALSE) {
  if (!is.null(color_col))
    check_metadata(ribo, color_col)
  
  coa_matrix <- extract_data(ribo, "count", position_to_rownames = TRUE,
                             only_annotated = only_annotated)
  
  coa_calculated <- .compute_coa(coa_matrix)
  
  if (object_only) {
    return(coa_calculated)
  }
  
  return(.plot_coa(coa_calculated, ribo[["metadata"]],
                   color_col, axes = axes, title, subtitle))
  
}

#' (internal) Correspondence Analysis computation function 
#'
#' @param raw_counts A matrix of counts, as exported by extract_data()
#' @return A coa object
#' @keywords internal
#' 
.compute_coa <- function(raw_counts = NULL) {
  res_coa <- ade4::dudi.coa(raw_counts[stats::complete.cases(raw_counts),
  ], scannf = FALSE, nf = 5)
  
  return(res_coa)
}


#' (internal) Plot a coa object 
#'
#' @param dudi.coa A coa object
#' @param metadata RiboClass's metadata dataframe
#' @inheritParams plot_coa
#' @return A ggplot object containing the COA
#' @keywords internal
#' 
.plot_coa <- function(dudi.coa = NULL,
                      metadata = NULL, color_col = NULL,
                      axes = c(1, 2), title = "default",
                      subtitle = "default") {
  if (is.null(color_col)) {
    color_colname <- color_col <- "Black"
  } else {
    color_colname <- color_col
    color_col <- as.factor(metadata[, color_col])
  }
  
  
  if (title == "default") {
    title <- "Correspondance analysis from count data"
  }
  
  if (subtitle == "default") {
    subtitle <- paste(ncol(dudi.coa$tab), "samples and", nrow(dudi.coa$tab),
                      "positions (all)")
  }
  
  plot.coa <- factoextra::fviz_ca_col(dudi.coa,
                                      repel = TRUE, col.col = color_col,
                                      pointsize = 2, labelsize = 4,
                                      axes = axes, title = title)
  
  plot.coa <- plot.coa + ggplot2::theme(text = ggplot2::element_text(size = 12))
  plot.coa <- plot.coa + ggplot2::labs(color = paste(color_colname), 
                                       subtitle = subtitle) +
    ggplot2::ylab(paste("Axis", axes[2], ":",
                        round(dudi.coa$eig[axes[2]]/sum(dudi.coa$eig) * 100, 1),
                        "%")) +
    ggplot2::xlab(paste("Axis", axes[1], ":",
                        round(dudi.coa$eig[axes[1]]/sum(dudi.coa$eig) * 100, 1),
                        "%"))
  return(plot.coa)
}
