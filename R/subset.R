
#' Keep only a subset of positions for all samples.
#'
#' @param ribo a riboClass object
#' @param positions_to_keep either a dataframe with a column "named_position" or a vector of named positions.
#'
#' @return a riboClass object with only selected positions for each sample.
#' @export
#'
#' @examples subset_ribo(ribo, c("18S_15","5.8S_4")) # returns a riboClass with only 2 selected positions
subset_ribo <- function(ribo, positions_to_keep) {
  #TODO : accept df without named_position + our df
  ribo_data <- ribo[["counts"]]
  position_df <- NULL
  if(is.data.frame(positions_to_keep)) {
    positions_df <- positions_to_keep
    positions_to_keep <- positions_to_keep[["named_position"]]
    
  }
  
  #check if positions_to_keep contains only existing positions. Throw a warning if not.
  existing_positions <- ribo[[1]][[1]][["named_position"]]
  if (!all(positions_to_keep %in% existing_positions)) {
    warning(paste(
      "One or more of the given positions to keep do not exist :",
      paste(positions_to_keep[!(positions_to_keep %in% existing_positions)], collapse = "; ")
    ))
  }
    
  #iterates through samples and keep only specified positions
  
  subsetted_ribo_data <- lapply(ribo_data, function(x) {
    if (!is.null(position_df)) {
      if ("renamed_position" %in% colnames(positions_df)) {
        x["renamed_position"] <-
          positions_df[["renamed_position"]][match(x[["named_position"]], positions_df[["named_position"]])]
        x["named_position"] <- x["renamed_position"]
        x <- subset(x, select = -c(renamed_position))
      }
      else {
        x <- x[which(x[,"named_position"] %in% positions_to_keep),]
      }
    }
    else {
      x <- x[which(x[,"named_position"] %in% positions_to_keep),]
      
    }
    
    return(x)
  })
  
  ribo[["counts"]] <- subsetted_ribo_data
  return(ribo)

  
}