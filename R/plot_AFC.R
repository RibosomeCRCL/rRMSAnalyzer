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
.calculate_afc <- function(raw.counts = NULL, metadata = NULL, order.by.col = NULL){
  ##### raw.counts = matrix with the counts
  ##### metadata = data.frame with sample annotation
  ##### order.by.col = index of metadata data.frame. It will order the raw.counts
  #TODO : better names
  res.coa <- ade4::dudi.coa(raw.counts[complete.cases(raw.counts), match(metadata[,order.by.col], colnames(raw.counts))], 
                            scannf = FALSE, 
                            nf = 5)
  
  return(res.coa)
  
}

#' Title
#'
#' @param ribo 
#' @param col_to_plot 
#' @param order_by_col 
#' @param col_for_color 
#'
#' @return
#' @export
#'
#' @examples
plot_AFC <- function(ribo,col_to_plot,order_by_col,col_for_color = NULL) {
  
  afc_matrix <- aggregate_samples_by_col(ribo[["raw_counts"]],col_to_plot,position_to_rownames = T)
   
  afc_calculated <- .calculate_afc(afc_matrix, ribo[["metadata"]],order_by_col)
  
  return(.plot_afc(afc_calculated,ribo[["metadata"]],col_for_color))
  
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
.plot_afc <- function(dudi.coa = NULL, metadata = NULL, col.by.col = NULL) {
  #### dudi.coa = output of dudi.coa function
  ##### metadata = data.frame with sample annotation
  ##### col.by.col = index of metadata data.frame. It will colour the raw.counts
  #TODO : not forcing metadata. 
  
  if(is.null(col.by.col)) {col.by.col = "Red"}
  else {col.by.col <- as.factor(metadata[,col.by.col])}
  
  plot.afc <- fviz_ca_col(dudi.coa, 
                          repel = T, 
                          col.col = col.by.col,
                          pointsize = 2,  #TODO : reduce size
                          labelsize = 4, 
                          axes = c(1,2),title = paste("Correspondence analysis of the raw counts on all genomic positions. Colored by: ", as.character(colnames(metadata)[col.by.col]))) + theme(text = element_text(size = 18))
  
  return(plot.afc) #TODO : find good colors and forms
  #TODO : 20 differents colors
  
} 