#' Plot a boxplot of a RiboClass object’s counts. 
#' @description This plot is useful to check if the samples are alike in their raw counts.
#' @param ribo A RiboClass object.
#' @param color_col Name of the column in the metadata used for coloring samples.
#' @param outlier Show boxplot outlier values.
#' 
#' @return A ggplot object.
#' @export
#'
#' @examples 
#' data("ribo_toy")
#' boxplot_count(ribo_toy,"run")
#' 
boxplot_count <- function(ribo, color_col=NA,outlier = T) {
  ribo_matrix <- extract_data(ribo,"count",position_to_rownames = T)
  return(.plot_boxplot_samples(ribo_matrix,"count","log10(count)",ribo[["metadata"]],color_col,outlier))

  
}

#' Title
#'
#' @param matrix 
#' @param values_col_name 
#' @param values_to_plot 
#' @param metadata 
#' @param color_col 
#' @param outlier 
#'
#' @return
#'
#' @examples
.plot_boxplot_sites <- function(matrix,values_col_name,values_to_plot, metadata,order_by,outlier, use_renamed_positions = T) {
  id_vars <- "sites"
  
  matrix_melted <- reshape2::melt(matrix, id.vars = id_vars,value.name = values_col_name,)
  
  shape_outlier <- NA
  if(outlier) shape_outlier <- 19
  p <- ggplot2::ggplot( matrix_melted, ggplot2::aes(x = reorder(named_position,!!sym(values_to_plot),na.rm = T), y = !!sym(values_to_plot))) +
    geom_boxplot(outlier.shape = shape_outlier) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  return(p)
}

#' (internal) plot boxplot for a given matrix of values
#'
#' @param matrix 
#' @param metadata 
#' @param color_col 
#' @param outlier 
#' @param values_col_name 
#' @param values_to_plot 
#'
#' @return
#'
#' @examples
.plot_boxplot_samples <- function(matrix,values_col_name,values_to_plot, metadata,color_col=NA,outlier) {
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

  return(p)
}
