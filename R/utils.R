check_metadata <- function(ribo,metadata_name) {
  
  if(!(metadata_name %in% colnames(ribo[["metadata"]]))) {
    stop("The metadata column specified (",metadata_name,") does not exist.\n Available columns : ",
         paste(colnames(ribo[["metadata"]]),collapse = ", "))
  }
}