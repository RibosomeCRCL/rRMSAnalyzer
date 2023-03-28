
#' Compute Relative Log Expression (RLE)
#'
#' @param mat Data matrix to compute RLE on
#'
#' @return a RLE matrix
#' @keywords internal
#' 
.compute_RLE <- function(mat = NULL){
  stopifnot(!is.null(mat))
  mat <- mat + 1
  med <- apply(mat, 1, stats::median)
  rle <- log2(mat/med)
  rle <- tidyr::gather(as.data.frame(rle))
  return(rle)
}


#' (internal) Compute and plot Relative Log Expression for a chosen data in as RiboClass
#'
#' @param mat matrix
#' @param show_outlier Show outliers on the plot ?
#' @param col_to_plot  data column to plot values from
#'
#' @return a ggplot object
#' @keywords internal
#' 
.plot_rle <- function(mat = NULL,show_outlier,col_to_plot) {
  
  # Non Standard Evaluation warning fix
  key <- NULL
  value <- NULL
  
  stopifnot(!is.null(mat))
  mat <- as.data.frame(mat)
  rle_calc <- .compute_RLE(mat = mat)
  outlier_shape=NA
  if(show_outlier) outlier_shape = 19
  mad <- - 2*mad(rle_calc$value, na.rm = T)
  rle_calc[["key"]] <- factor(rle_calc[["key"]],levels = unique(rle_calc[["key"]]))
  rle_grouped <- rle_calc %>% dplyr::group_by(key) %>% dplyr::summarise(median = stats::median(value, na.rm = TRUE))
  
  rle.plot <- ggplot2::ggplot(rle_calc, ggplot2::aes(x = key, value)) + 
    ggplot2::geom_boxplot(outlier.shape = outlier_shape, outlier.size = 0.3,
                 fill = ifelse(rle_grouped$median < mad, "red","white")) +
    ggplot2::ylab(paste0("RLE")) +
    ggplot2::xlab("Sample") +
    ggplot2::ggtitle("RLE plot")+
    ggplot2::geom_hline(yintercept = 2*mad(rle_calc$value, na.rm = T), colour = "blue") + 
    ggplot2::geom_hline(yintercept = - 2*mad(rle_calc$value, na.rm = T), colour = "blue") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))
 
   return(rle.plot)
}

#' Relative Log Expression plot of a RiboClass object’s counts. 
#'
#' @param ribo a RiboClass object.
#' @param outlier Show boxplot outlier values. 
#'
#' @return A ggplot object. Samples with a median lower than median(RLE counts)-2*MAD (Median Absolute Deviation) are colored in red.
#' @export
#'
#' @examples 
#' data("ribo_toy")
#' plot_rle(ribo_toy)
plot_rle <- function(ribo, show_outlier=TRUE) {
  
  rle_matrix <- extract_data(ribo, col = "count",position_to_rownames = T)
  return(.plot_rle(rle_matrix,show_outlier = show_outlier,"count"))
}
