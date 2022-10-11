#' Rename RNAs in your RiboClass
#'
#'  
#' 
#' @param ribo a RiboClass object
#' @param new_names the new names for your RNA (by order of rna size)
#'
#' 
#'
#' @return a RiboClass with updated rna names in data.
#' @export
#'
#' @examples
#' data("ribo_toy")
#' ribo_toy <- rename_rna(ribo_toy ,c("5S","5.8S","18S","28S"))

rename_rna <- function(ribo,new_names=c("5S","5.8S","18S","28S")) {
  
  sample_list <- ribo[["data"]]
  rna_names <- ribo[["rna_names"]]
  
  if(nrow(rna_names) != length(new_names)) {
    stop("Different numbers of RNA names in your RiboClass (",nrow(rna_names),") and the list given (",length(new_names),").")
  }
  rna_names[3] <- new_names
  
  #change the RNA names inside each sample
  sample_list_renamed <- lapply(sample_list, function(x) {
    x[,1] <- rna_names[,3][match(x[,1], rna_names[,2])]
    return(x)
  })
  rna_names[2] <- rna_names[3]
  rna_names <- rna_names[,1:2]
  # Update nomenclature according to the new RNAs
  sample_list_renamed <- .generate_riboclass_named_position(sample_list_renamed,1,2)
  
  ribo[["data"]] <- sample_list_renamed
  ribo[["rna_names"]] <- rna_names
  return(ribo)
}