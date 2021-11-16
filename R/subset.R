
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
  
  ribo_data <- ribo[["counts"]]
  
  #check if positions_to_keep contains only existing positions. Throw a warning if not.
  
  existing_positions <- ribo[[1]][[1]][["named_position"]]
  if (!isTRUE(positions_to_keep %in% existing_positions)) {
    warning(paste(
      "One or more of the given positions to keep do not exist :",
      paste(positions_to_keep[!(positions_to_keep %in% existing_positions)], collapse = "; ")
    ))
  }
    
  #iterates through samples and keep only specified positions
  
  subsetted_ribo_data <- lapply(ribo_data, function(x) {
    x <- x[which(x[,"named_position"] %in% positions_to_keep),]
    return(x)
  })
  
  ribo[["counts"]] <- subsetted_ribo_data
  return(ribo)

  
}