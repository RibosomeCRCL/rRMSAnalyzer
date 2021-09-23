
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
  library(tidyr)
  rle <- gather(as.data.frame(rle))
  return(rle)
}


#' Title
#'
#' @param mat 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
.plot_RLE <- function(mat = NULL, ...) {
  stopifnot(!is.null(mat))
  mat <- as.data.frame(mat)
  rle.calc <- .calculate_relative_log_expression(mat = mat)
  library(ggplot2)
  rle.plot <- ggplot(rle.calc, aes(x = key, value)) + geom_boxplot(outlier.size = 0.3, fill = "red") + ylab("log2(C-score/median)") + xlab("Sample ID") + ggtitle("RLE plot") + theme(axis.text.x = element_text(angle=45, vjust = 0.8)) +  geom_hline(yintercept = 2*mad(rle.calc$value, na.rm = T), colour = "blue") + geom_hline(yintercept = - 2*mad(rle.calc$value, na.rm = T), colour = "blue")
  return(rle.plot)
}

#' Title
#'
#' @param ribo 
#' @param col_to_plot 
#'
#' @return
#' @export
#'
#' @examples
plot_RLE <- function(ribo, col_to_plot) {
  
  rle_matrix <- aggregate_samples_by_col(ribo[["raw_counts"]],col_to_keep = col_to_plot,position_to_rownames = T)
  return(.plot_RLE(rle_matrix))
}
