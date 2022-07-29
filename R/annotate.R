#' Annotate sites according to a dataframe containing annotation
#' 
#' This function will fill the site column with a nomenclature given in anodf. 
#' 
#' 
#' @details Your dataframe must contain
#' \itemize{
#'  \item{siteID : the column containing the nomenclature for each site e.g "Am26"}
#'  \item{rna : name of the RNA where the site is located e.g "18S"}
#'  \item{position : site position on the RNA e.g 234} 
#' }
#' Other columns will be ignored. 

#' 
#' @param ribo the riboclass in which you want to annotate sites
#' @param anodf the dataframe containing your annotations
#'
#' @return
#' @export
#'
#' @examples
annotate_site <- function(ribo, anodf) {
  ribo_data <- ribo[["data"]]
  
  
  #check if anodf contains only existing positions. Throw a warning if not.
  existing_positions <- ribo[[1]][[1]][["named_position"]]
  if (!all(anodf[["named_position"]] %in% existing_positions)) {
    warning(paste(
      "One or more of the given positions to keep do not exist :",
      paste(anodf["named_position"][!(anodf["named_position"] %in% existing_positions)], collapse = "; ")
    ))
  }
  
  if("site" %in% colnames(anodf)) {
  
  subsetted_ribo_data <- lapply(ribo_data, function(x) {
    x["site"] <-  anodf[["siteID"]][match(x[["named_position"]], anodf[["named_position"]])]
    
    
    return(x)
  })
  }
  else {
    stop("the \"site\" column is missing in annotation data. Please add it with site names.")
  }
  
  ribo[["data"]] <- subsetted_ribo_data
  return(ribo)
  
}