#' Annotate sites according to a dataframe containing annotation
#' 
#' Mark sites of interest by giving them a custom name. 
#' Some analyses and plots can be applied on these specific sites when only_annotated parameter 
#' is available.
#' 
#' @details
#' This function will fill the "site" column in your riboclass's data with a nomenclature given in annot. 
#' 
#' @param ribo a RiboClass object to annotate, see : 
#' \code{\link{load_ribodata}}
#' @param annot The dataframe containing annotations.
#' @param anno_rna Name or index of the column in annot containing RNAs' name.
#' @param anno_pos Name or index of the column in annot containing site position inside RNA.
#' @param anno_value Name or index of the column in annot containing nomenclature to apply.
#' 
#' @return An annotated RiboClass (the site column in data should be filled with names from annot_value for known position).
#' @export
#'
#' @examples
#' data("ribo_toy")
#' data("human_methylated")
#' ribo_toy <- rename_rna(ribo_toy ,c("5S","5.8S","18S","28S"))
#' ribo_toy <- annotate_site(ribo_toy,human_methylated,anno_value ="Nomenclature")
#'
annotate_site <- function(ribo, annot, anno_rna = 2, anno_pos = 1, anno_value = 3) {
  ribo_data <- ribo[["data"]]
  
  if(!("named_position" %in% colnames(annot))) {
    annot <- .generate_name_positions(annot,anno_rna,anno_pos)
  }
  
  #Check if annot has the RNA as the RiboClass
 if(sum(!(annot[[anno_rna]] %in% ribo[["rna_names"]][["current_name"]])) != 0){
   stop("Mismatch between annot's RNA names and RiboClass's RNA names\n you can rename your RiboClass' RNAs with rename_rna")
 }
  
  #check if annot contains only existing positions. Throw a warning if not.
  existing_positions <- ribo[[1]][[1]][["named_position"]]
  if (!all(annot[["named_position"]] %in% existing_positions)) {
    warning(paste(
      "One or more of the given positions to annotate do not exist :",
      paste(annot[anno_value][!(annot["named_position"] %in% existing_positions)], collapse = "; ")
    ))
  }
  
  if("site" %in% colnames(annot)) {
    anno_value <- "site"
  }
  
  else if (is.null(anno_value)){
    stop("the \"site\" column is missing in annotation data. Please add it with site names.")
  }

  subsetted_ribo_data <- lapply(ribo_data, function(x) {
    x["site"] <-  annot[[anno_value]][match(x[["named_position"]], annot[["named_position"]])]
    return(x)
  })
  
  ribo[["data"]] <- subsetted_ribo_data
  return(ribo)
  
}