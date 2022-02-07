#' (internal) factorial correspondence analysis calculation function 
#'
#' @param raw.counts 
#' @param metadata 
#' @param order.by.col 
#'
#' @return a coa object
#' @export
#'
#' @examples .calculate_fca(raw_counts = cscore_matrix)
.calculate_fca <- function(raw.counts = NULL){
  #TODO : better names
  res.coa <- ade4::dudi.coa(raw.counts[complete.cases(raw.counts),], 
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
#' @param axis COA axis to plot. First and second axis will be plotted by default
#'
#' @return a FCA plot
#' @export
#'
#' @examples plot_fca(ribo = ribo, col_for_color = "condition")
plot_fca <- function(ribo,col_for_color = NULL, axis = c(1,2)) {
  
  fca_matrix <- aggregate_samples_by_col(ribo[["counts"]],"Count",position_to_rownames = T)
   
  fca_calculated <- .calculate_fca(fca_matrix)
  
  return(.plot_fca(fca_calculated,ribo[["metadata"]],col_for_color, axis = axis))
  
}

#' (internal) Plot a coa object 
#'
#' @param dudi.coa a coa object
#' @param metadata 
#' @param col.by.col 
#'
#' @return
#' @export
#'
#' @examples
.plot_fca <- function(dudi.coa = NULL, metadata = NULL, col.by.col = NULL, axis = c(1,2)) {
  #TODO : not forcing metadata. 
  if(is.null(col.by.col)) {col.by.col = "Red"}
  else {col.by.col <- as.factor(metadata[,col.by.col])}
  
  plot.fca <- fviz_ca_col(dudi.coa, 
                          repel = T,
                          col.col = col.by.col,
                          pointsize = 2,  #TODO : reduce size
                          labelsize = 4, 
                          axes = axis,
                          title = paste("Correspondence analysis of the raw counts on all genomic positions (", nrow(fca$tab),")")) 
  plot.fca <- plot.fca + theme(text = element_text(size = 12)) 
  plot.fca <- plot.fca + labs(color = paste(col.by.col), subtitle = paste(ncol(fca$tab), "samples"))  # col.by.col returns "e"

  
  return(plot.fca) #TODO : find good colors and forms
  #TODO : 20 differents colors
  
} 
