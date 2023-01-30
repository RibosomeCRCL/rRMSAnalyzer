#' Update count values inside a RiboClass with a matrix of values
#' 
#' @param ribo a RiboClass object
#' @param update_matrix a position x sample matrix containing the new values 
#'
#' @return a RiboClass with updated values
#' @keywords internal
#'
.update_ribo_count_with_matrix <- function(ribo, update_matrix) {
  #first, check if we have the sample name in our column
  update_df <- as.data.frame(update_matrix)
  col_names <- sort(names(update_df))
  riboclass_names <- sort(names(ribo[["data"]]))
  
  if(!identical(col_names,riboclass_names)) {
    stop("mismatch between samplenames and matrix's samples names")
  }
  
  # For each sample in the RiboClass, replace count values with matrix's ones
  
  count_list <- ribo[["data"]]
  
  for(sample in names(count_list)) {
    count_list[[sample]]["count"] <- update_df[sample]
  }
  
  ribo[["data"]] <- count_list
  
  return(ribo)
  
  
}