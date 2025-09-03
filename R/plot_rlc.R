
#' Compute Relative Log Expression (RLC)
#'
#' @param mat Data matrix to compute RLC on
#'
#' @return a RLC matrix
#' @keywords internal
#' 
.compute_RLC <- function(mat = NULL){
  stopifnot(!is.null(mat))
  mat <- mat + 1
  med <- apply(mat, 1, stats::median)
  rlc <- log2(mat/med)
  rlc <- tidyr::gather(as.data.frame(rlc))
  return(rlc)
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
.plot_rlc <- function(mat = NULL, show_outlier, col_to_plot) {
  
  # Non Standard Evaluation warning fix
  key <- NULL
  value <- NULL
  
  stopifnot(!is.null(mat))
  mat <- as.data.frame(mat)
  rlc_calc <- .compute_RLC(mat = mat)
  outlier_shape <- NA
  if(show_outlier) outlier_shape <- 19
  mad <- - 2*mad(rlc_calc$value, na.rm = TRUE)
  rlc_calc[["key"]] <- factor(rlc_calc[["key"]],levels = unique(rlc_calc[["key"]]))
  rlc_grouped <- rlc_calc %>% dplyr::group_by(key) %>% dplyr::summarise(median = stats::median(value, na.rm = TRUE))
  
  rlc.plot <- ggplot2::ggplot(rlc_calc, ggplot2::aes(x = key, value)) + 
    ggplot2::geom_boxplot(outlier.shape = outlier_shape, outlier.size = 0.3,
                 fill = ifelse(rlc_grouped$median < mad, "red","white")) +
    ggplot2::ylab(paste0("RLC")) +
    ggplot2::xlab("Sample") +
    ggplot2::ggtitle("RLC plot")+
    ggplot2::geom_hline(yintercept = 2*mad(rlc_calc$value, na.rm = TRUE),
                        colour = "blue") + 
    ggplot2::geom_hline(yintercept = - 2*mad(rlc_calc$value, na.rm = TRUE),
                        colour = "blue") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))
 
   return(rlc.plot)
}

#' Relative Log Expression plot of a RiboClass object’s counts. 
#'
#' @param ribo a RiboClass object.
#' @param show_outlier Show boxplot outlier values. 
#'
#' @return A ggplot object. Samples with a median lower than median(RLC counts)-2*MAD (Median Absolute Deviation) are colored in red.
#' @export
#'
#' @examples 
#' data("ribo_toy")
#' plot_rlc(ribo_toy)
plot_rlc <- function(ribo, show_outlier=TRUE) {
  
  rlc_matrix <- extract_data(ribo, col = "count",position_to_rownames = TRUE)
  return(.plot_rlc(rlc_matrix,show_outlier = show_outlier,"count"))
}
