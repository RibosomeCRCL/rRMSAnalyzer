#' Title
#'
#' @param raw.counts 
#' @param metadata 
#' @param order.by.col 
#'
#' @return
#' @export
#'
#' @examples
.calculate_fca <- function(raw.counts = NULL, metadata = NULL, order.by.col = NULL){
  ##### raw.counts = matrix with the counts
  ##### metadata = data.frame with sample annotation
  ##### order.by.col = index of metadata data.frame. It will order the raw.counts
  #TODO : better names
  res.coa <- ade4::dudi.coa(raw.counts[complete.cases(raw.counts), match(metadata[,order.by.col], colnames(raw.counts))], 
                            scannf = FALSE, 
                            nf = 5)
  
  return(res.coa)
  
}

#' factorial correspondence analysis of a riboclass for a given column
#'
#' @param ribo A RiboClass object
#' @param col_to_plot column containing the values to plot the fca on
#' @param order_by_col column containing the sample names for linking data and metadata
#' @param col_for_color metadata column to colorize samples
#'
#' @return a FCA plot
#' @export
#'
#' @examples
plot_fca <- function(ribo,col_to_plot,order_by_col,col_for_color = NULL) {
  
  fca_matrix <- aggregate_samples_by_col(ribo[["raw_counts"]],col_to_plot,position_to_rownames = T)
   
  fca_calculated <- .calculate_fca(fca_matrix, ribo[["metadata"]],order_by_col)
  
  return(.plot_fca(fca_calculated,ribo[["metadata"]],col_for_color))
  
}

#' Title
#'
#' @param dudi.coa 
#' @param metadata 
#' @param col.by.col 
#'
#' @return
#' @export
#'
#' @examples
.plot_fca <- function(dudi.coa = NULL, metadata = NULL, col.by.col = NULL) {
  #### dudi.coa = output of dudi.coa function
  ##### metadata = data.frame with sample annotation
  ##### col.by.col = index of metadata data.frame. It will colour the raw.counts
  #TODO : not forcing metadata. 
  if(is.null(col.by.col)) {col.by.col = "Red"}
  else {col.by.col <- as.factor(metadata[,col.by.col])}
  
  plot.fca <- fviz_ca_col(dudi.coa, 
                          repel = T,
                          col.col = col.by.col,
                          pointsize = 2,  #TODO : reduce size
                          labelsize = 4, 
                          axes = c(1,2),title = paste("Correspondence analysis of the raw counts on all genomic positions")) + theme(text = element_text(size = 12)) + labs(color = paste(col.by.col), subtitle = paste(length(metadata[,"filename"]), "samples"))  # col.by.col returns "e"

  
  return(plot.fca) #TODO : find good colors and forms
  #TODO : 20 differents colors
  
} 