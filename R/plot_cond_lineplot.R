#' Title
#'
#' @param ribo 
#' @param value_col 
#' @param metadata_col 
#'
#' @return
#' @export
#'
#' @examples
plot_lineplot <- function(ribo,value_col,metadata_col) {
 ribo_concat <- mean_samples_by_conditon(ribo,value_col,metadata_col )
  
 lineplot <- ggplot(ribo_concat, aes(x=named_position, y=mean, colour=!!sym(metadata_col), group = !!sym(metadata_col))) + 
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1, position = position_dodge(0.1)) +
    geom_line(position = position_dodge(0.1)) +
    geom_point(position = position_dodge(0.1), size=3, shape=21, fill="white") + 
    rotate_x_text(angle = 90)
 
 return(lineplot)
}