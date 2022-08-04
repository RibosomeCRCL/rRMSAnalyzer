#' annotate sites according to a dataframe containing annotation
#' 
#' This function will fill the "site" column in your riboclass's data with a nomenclature given in annodf. 
#' 
#' @param ribo the riboclass in which you want to annoate sites
#' @param annodf the dataframe containing your annoations
#' @param anno_rna name or position of the column in annodf containing RNAs' name
#' @param anno_pos name or position of the column in annodf containing site position inside RNA.
#' @param anno_col name or position of the column in annodf containing nomenclature to apply
#' 
#' @return
#' @export
#'
#' @examples
#' ribo_toy <- rename_rna(ribo_toy ,c("5S","5.8S","18S","28S"))
#' ribo_toy <- annotate_site(ribo_toy,human_methylated,anno_col ="Nomenclature")
#'
annotate_site <- function(ribo, annodf, anno_rna = 2, anno_pos = 1 , anno_col = NULL) {
  ribo_data <- ribo[["data"]]
  
  if(!("named_position" %in% colnames(annodf))) {
    annodf <- generate_name_positions(annodf,anno_rna,anno_pos)
  }
  
  #check if annodf contains only existing positions. Throw a warning if not.
  existing_positions <- ribo[[1]][[1]][["named_position"]]
  if (!all(annodf[["named_position"]] %in% existing_positions)) {
    warning(paste(
      "One or more of the given positions to keep do not exist :",
      paste(annodf["named_position"][!(annodf["named_position"] %in% existing_positions)], collapse = "; ")
    ))
  }
  
  if("site" %in% colnames(annodf)) {
    anno_col <- "site"
  }
  
  else if (is.null(anno_col)){
    stop("the \"site\" column is missing in annotation data. Please add it with site names.")
    
  }

  subsetted_ribo_data <- lapply(ribo_data, function(x) {
    x["site"] <-  annodf[[anno_col]][match(x[["named_position"]], annodf[["named_position"]])]
    
    return(x)
  })
  
  ribo[["data"]] <- subsetted_ribo_data
  return(ribo)
  
}