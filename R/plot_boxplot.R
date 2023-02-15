#' Plot a boxplot of a RiboClass object’s counts. 
#' @description This plot is useful to check if the samples are alike in their raw counts.
#' @param ribo A RiboClass object.
#' @param color_col Name of the column in the metadata used for coloring samples.
#' @param outlier Show boxplot outlier values.
#' @param horizontal Show boxplot horizontally.
#' @return A ggplot object.
#' @export
#'
#' @examples 
#' data("ribo_toy")
#' boxplot_count(ribo_toy,"run")
#' 
boxplot_count <- function(ribo, color_col=NA,outlier = TRUE, horizontal = FALSE) {
  ribo_matrix <- extract_data(ribo,"count",position_to_rownames = TRUE)
  return(.plot_boxplot_samples(ribo_matrix,"count","log10(count)",
                               ribo[["metadata"]],color_col,outlier,horizontal = horizontal))

  
}

#' Plot boxplot of annotated sites c-score.
#' Sites are sorted by 
#' @inheritParams boxplot_count 
#'
#' @return
#' @export
#'
#' @examples
boxplot_cscores <- function(ribo,outlier = TRUE, horizontal = FALSE) {
  
  ribo_m <- extract_data(ribo,only_annotated = TRUE)
  only_annotated = TRUE
  
  if(nrow(ribo_m) == 0) {
    stop("No annotated site found. Please use annotate_site() on your RiboClass before calling this function.")
    }
  
  return(.plot_boxplot_sites(ribo_m,
                      values_to_plot = "cscore",outlier = TRUE,
                      horizontal = horizontal))
}

#' Title
#'
#' @param matrix 
#' @param values_col_name 
#' @param values_to_plot 
#' @param outlier 
#'
#' @return
#'
.plot_boxplot_sites <- function(matrix,values_to_plot,outlier, horizontal) {
  id_vars <- "site"

  
  matrix_melted <- reshape2::melt(matrix, id.vars = id_vars,
                                  value.name = values_to_plot)
  
  shape_outlier <- NA
  if(outlier) shape_outlier <- 19
  p <- ggplot2::ggplot( matrix_melted, ggplot2::aes(
    x = reorder(site,!!rlang::sym(values_to_plot),na.rm = T),
    y = !!rlang::sym(values_to_plot))) +
    ggplot2::geom_boxplot(outlier.shape = shape_outlier) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggplot2::xlab("Site") + ggplot2::ylab("C-score")
  
  if(horizontal) p <- p + ggplot2::coord_flip()
  
  return(p)
}

#' (internal) plot boxplot for a given matrix of values
#'
#' @param matrix Matrix of values
#' @param metadata 
#' @param color_col 
#' @param outlier 
#' @param values_col_name 
#' @param values_to_plot 
#'
#' @return ggplot boxplot
#'
.plot_boxplot_samples <- function(matrix,values_col_name,values_to_plot,
                                  metadata,color_col=NA,outlier, horizontal) {
  id_vars <- "Sample"
  matrix_inv <- as.data.frame(t(matrix))
  matrix_inv <- tibble::rownames_to_column(matrix_inv, "Sample")
  
  if(!is.na(color_col)){
    matrix_inv <- cbind(matrix_inv,metadata[color_col])
    id_vars <- c(color_col,"Sample")
    } #TODO : proper merge
  
  matrix_melted <- reshape2::melt(matrix_inv, id.vars = id_vars,value.name = values_col_name)
  

  shape_outlier <- NA
  if(outlier) shape_outlier <- 19
  
  if(is.na(color_col)) {
    p <- ggplot2::ggplot(matrix_melted, ggplot2::aes_string(x = "Sample", y = values_to_plot)) +
      ggplot2::geom_boxplot(outlier.shape = shape_outlier) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1))
   
  }
  else {
  p <- ggplot2::ggplot(matrix_melted, ggplot2::aes_string(x = "Sample", y = values_to_plot,fill = color_col)) +
    ggplot2::geom_boxplot(outlier.shape = shape_outlier) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1))
 
  }
 p <-  p +  ggplot2::geom_hline(yintercept = 2,colour = "blue") 
  if(horizontal) p <- p + ggplot2::coord_flip()
 
  return(p)
}
