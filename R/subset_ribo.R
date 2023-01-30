#' Keep only a subset of positions for all samples.
#'
#' @param ribo a RiboClass object, see constructor : 
#' \code{\link{create_riboclass}}
#' @param positions_to_keep either a dataframe with a column "named_position" or a vector of named positions.
#' @param anno_rna name or index of the column containing RNA names in annotation data
#' @param anno_pos name or index of the column containing position in annotation data
#'
#' @return a RiboClass object with only selected positions for each sample.
#' @export
#'
#' @examples 
#' data("ribo_toy")
#' subset_ribo(ribo_toy, c("NR_046235.3_18S_0681","NR_046235.3_18S_0682"))
subset_ribo <- function(ribo, positions_to_keep, anno_rna, anno_pos) {
  #TODO : accept df without named_position + our df
  ribo_data <- ribo[["data"]]
  positions_df <- NULL
  if(is.data.frame(positions_to_keep)) {
    if(!("named_position" %in% positions_to_keep)) { 
      positions_to_keep <- .generate_name_positions(positions_to_keep,anno_rna,anno_pos)
    }
    
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
    # fix R CMD check note due to the non-standard evaluation in the subset function call
    renamed_position <- NULL 
    if (!is.null(positions_df)) {
      x <- x[which(x[,"named_position"] %in% positions_df[,"named_position"]),]
      
      if ("renamed_position" %in% colnames(positions_df)) {
        x["renamed_position"] <-  positions_df[["renamed_position"]][match(x[["named_position"]], positions_df[["named_position"]])]
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
  
  ribo[["data"]] <- subsetted_ribo_data
  return(ribo)

  
}