#' Plot a boxplot of all counts for each sample
#' @description This plot is useful to check if the samples are alike in their raw counts.
#' @param ribo a riboclass object containing your samples
#' @param metadata_col (optional) metadata column to colorize your boxes
#'
#' @return a ggplot object
#' @export
#'
#' @examples plot_boxplot_count(ribo,"group")
plot_boxplot_count <- function(ribo, metadata_col=NA,show_outlier = T) {
  ribo_matrix <- aggregate_samples_by_col(ribo[[1]],"Count",position_to_rownames = T)
  return(.plot_boxplot_samples(ribo_matrix,"count","log10(count)",ribo[["metadata"]],metadata_col,show_outlier))

  
}

#' Title
#'
#' @param matrix 
#' @param values_col_name 
#' @param values_to_plot 
#' @param metadata 
#' @param metadata_col 
#' @param show_outlier 
#'
#' @return
#' @export
#'
#' @examples
.plot_boxplot_sites <- function(matrix,values_col_name,values_to_plot, metadata,order_by,show_outlier, use_renamed_positions = T) {
  id_vars <- "named_position"
  
  matrix_melted <- reshape2::melt(matrix, id.vars = id_vars,value.name = values_col_name,)
  
  shape_outlier <- NA
  if(show_outlier) shape_outlier <- 19
  p <- ggplot2::ggplot( matrix_melted, aes(x = reorder(named_position,!!sym(values_to_plot),na.rm = T), y = !!sym(values_to_plot))) +
    geom_boxplot(outlier.shape = shape_outlier) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  return(p)
}

#' (internal) plot boxplot for a given matrix of values
#'
#' @param matrix 
#' @param metadata 
#' @param metadata_col 
#' @param show_outlier 
#' @param values_col_name 
#' @param values_to_plot 
#'
#' @return
#' @export
#'
#' @examples
.plot_boxplot_samples <- function(matrix,values_col_name,values_to_plot, metadata,metadata_col=NA,show_outlier) {
  id_vars <- "sample"
  matrix_inv <- as.data.frame(t(matrix))
  matrix_inv <- tibble::rownames_to_column(matrix_inv,"sample")
  
  if(!is.na(metadata_col)){
    matrix_inv <- cbind(matrix_inv,metadata[metadata_col])
    id_vars <- c(metadata_col,"sample")
    } #TODO : proper merge
  
  matrix_melted <- reshape2::melt(matrix_inv, id.vars = id_vars,value.name = values_col_name)
  

  shape_outlier <- NA
  if(show_outlier) shape_outlier <- 19
  
  if(is.na(metadata_col)) {
    p <- ggplot2::ggplot(matrix_melted, aes_string(x = "sample", y = values_to_plot)) +
      geom_boxplot(outlier.shape = shape_outlier) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
   
  }
  else {
  p <- ggplot2::ggplot(matrix_melted, aes_string(x = "sample", y = values_to_plot,fill = metadata_col)) +
    geom_boxplot(outlier.shape = shape_outlier) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
 
  }
 p <-  p +  geom_hline(yintercept = 2,colour = "blue") 

  return(p)
}