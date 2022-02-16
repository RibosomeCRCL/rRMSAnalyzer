
#' Title
#'
#' @param mat 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
.calculate_relative_log_expression <- function(mat = NULL, ...){
  stopifnot(!is.null(mat))
  mat <- mat + 1
  med <- apply(mat, 1, median)
  rle <- log2(mat/med)
  rle <- tidyr::gather(as.data.frame(rle))
  return(rle)
}


#' Title
#'
#' @param mat 
#' @param ... 
#' @param show_outlier 
#'
#' @return
#' @export
#'
#' @examples
.plot_RLE <- function(mat = NULL,show_outlier,col_to_plot) {
  stopifnot(!is.null(mat))
  mat <- as.data.frame(mat)
  rle.calc <- .calculate_relative_log_expression(mat = mat)
  outlier_shape=NA
  if(show_outlier) outlier_shape = 19
  mad <- - 2*mad(rle.calc$value, na.rm = T)
  rle_grouped <- rle.calc %>% group_by(key) %>% summarise(median = median(value, na.rm = TRUE))
  
  rle.plot <- ggplot(rle.calc, aes(x = key, value)) + 
    geom_boxplot(outlier.shape =outlier_shape ,outlier.size = 0.3, fill = ifelse(rle_grouped$median < mad, "red","white")) +
    ylab(paste0("log2(",col_to_plot,"/median)")) +
    xlab("Sample") +
    ggtitle("RLE plot") + 
    geom_hline(yintercept = 2*mad(rle.calc$value, na.rm = T), colour = "blue") + 
    geom_hline(yintercept = - 2*mad(rle.calc$value, na.rm = T), colour = "blue") +
    theme(axis.text.x = element_text(colour = ifelse(rle_grouped$median < mad,"red","black"),angle=90))
 
   return(rle.plot)
}

#' Title
#'
#' @param ribo 
#' @param col_to_plot 
#' @param show_outlier 
#'
#' @return
#' @export
#'
#' @examples
plot_RLE <- function(ribo, col_to_plot,show_outlier) {
  
  rle_matrix <- aggregate_samples_by_col(ribo[["counts"]],col_to_keep = col_to_plot,position_to_rownames = T)
  return(.plot_RLE(rle_matrix,show_outlier = show_outlier,col_to_plot))
}
