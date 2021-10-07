
#' Title
#'
#' @param ribo 
#' @param positions_to_keep 
#'
#' @return
#' @export
#'
#' @examples
subset_ribo <- function(ribo, positions_to_keep) {
  
  ribo_data <- ribo[["raw_counts"]]
  
  #iterates through samples and keep only specified positions
  
  subsetted_ribo_data <- lapply(ribo_data, function(x) {
    x <- x[which(x[,"named_position"] %in% positions_to_keep),]
    return(x)
  })
  
  ribo[["raw_counts"]] <- subsetted_ribo_data
  return(ribo)

  
}