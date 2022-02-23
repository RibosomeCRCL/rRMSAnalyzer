#' Title
#'
#' @param ribo 
#' @param values_col 
#' @param metadata_col 
#'
#' @return
#' @export
#'
#' @examples
plot_barpos <- function(ribo,values_col,metadata_col) {
  ribo_concat <- mean_samples_by_conditon(ribo,values_col,metadata_col)
  barpos <- ggplot(ribo_concat,  aes(x=!!sym(metadata_col), y=mean, fill = !!sym(metadata_col))) + geom_bar(stat="identity", alpha=0.8) + geom_errorbar( aes(x=!!sym(metadata_col), ymin=mean-sd, ymax=mean+sd)) + facet_wrap(~named_position) + rotate_x_text(angle = 90)
  return(barpos)
}