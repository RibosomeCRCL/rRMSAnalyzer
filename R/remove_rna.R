#' Remove one or more RNA among your RiboClass' samples
#'
#' @param ribo a RiboClass
#' @param rna_to_remove names of the RNA to remove
#'
#' @return
#' @export
#'
#' @examples
remove_rna <- function(ribo,rna_to_remove) {
  
  if(!(rna_to_remove %in% ribo[["rna_names"]][["current_name"]])) {
    stop("the RNA names given do not exist in the RiboClass")
  }
  
  ribo_data <- lapply(ribo[["data"]], function(x) {
    x <- x[which(x["rna"] != rna_to_remove),]
    x["rna"] <- droplevels(x["rna"]) # drop unused factor level
    return(x)
  })
  ribo[["data"]] <- ribo_data
  ribo[["rna_names"]] <- ribo[["rna_names"]][which(ribo[["rna_names"]]["current_name"] != rna_to_remove),]
  
  return(ribo)
}