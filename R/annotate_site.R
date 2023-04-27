#' Annotate sites according to a dataframe containing annotation
#' 
#' Annotate sites of interest by giving them a custom name. 
#' Some analyses and plots can be applied on these specific sites when only_annotated parameter 
#' is available.
#' 
#' @details
#' This function will fill the 'site' column in your RiboClass's data with a nomenclature given in annot. 
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
#' data('ribo_toy')
#' data('human_methylated')
#' ribo_toy <- rename_rna(ribo_toy ,c('5S','5.8S','18S','28S'))
#' ribo_toy <- annotate_site(ribo_toy,human_methylated,anno_value ='Nomenclature')
#'
annotate_site <- function(ribo, annot, anno_rna = 2, anno_pos = 1, anno_value = 3) {
  ribo_data <- ribo[["data"]]
  
  if (!("named_position" %in% colnames(annot))) {
    annot <- .generate_name_positions(annot, anno_rna, anno_pos)
  }
  
  # Check if annot has the RNA as the RiboClass
  if (sum(!(annot[[anno_rna]] %in% ribo[["rna_names"]][["current_name"]])) != 0) {
    stop("Mismatch between annot's RNA names and RiboClass's RNA names. You can rename your RiboClass' RNAs with rename_rna")
  }
  
  # check if annot contains only existing positions. Throw a
  # warning if not.
  existing_positions <- ribo[[1]][[1]][["named_position"]]
  if (!all(annot[["named_position"]] %in% existing_positions)) {
    warning(paste("One or more of the given positions to annotate do not exist :",
                  paste(annot[anno_value][!(annot["named_position"] %in% existing_positions)],
                        collapse = "; ")))
  }
  
  if ("site" %in% colnames(annot)) {
    anno_value <- "site"
  } else if (is.null(anno_value)) {
    stop("the \"site\" column is missing in annotation data. Please add it with site names.")
  }
  
  subsetted_ribo_data <- lapply(ribo_data, function(x) {
    x["site"] <- annot[[anno_value]][match(x[["named_position"]], annot[["named_position"]])]
    return(x)
  })
  
  ribo[["data"]] <- subsetted_ribo_data
  return(ribo)
  
}

#' Remove site annotations of a given RiboClass
#'
#' @param ribo A RiboClass object.
#'
#' @return A RiboClass object where all annotated position have been replaced by
#' NA, the default value.
#' @export
#'
#' @examples
#' data('ribo_toy')
#' remove_annotation(ribo_toy)
remove_annotation <- function(ribo) {
  ribo[["data"]] <- lapply(ribo[["data"]], function(x) {
    x["site"] <- NA
    return(x)
  })
  
  return(ribo)
}

#' Keep only a subset of the current annotation
#'
#' @param ribo a RiboClass
#' @param annotation_to_keep vector containing annotated sites'name to keep
#'
#' @return a RiboClass where only annotated sites within annotation_to_keep are
#' still annotated.
#' @export
#'
#' @examples
#' data('ribo_toy')
#' data('human_methylated')
#' ribo_toy <- rename_rna(ribo_toy)
#' ribo_toy <- annotate_site(ribo_toy,human_methylated)
#' ribo_toy <- keep_selected_annotation(ribo_toy, c("28S_Am1310","28S_Cm2848"))
#' 
keep_selected_annotation <- function(ribo, annotation_to_keep) {
  
  sample_1 <- ribo[["data"]][[1]]
  
  current_annotation <- sample_1[which(!is.na(sample_1[["site"]])),
                                 c("rnapos","rna","site")]
  
  ribo <- remove_annotation(ribo)
  
  new_annotation <- current_annotation[
    which(current_annotation[["site"]] %in% annotation_to_keep),]
  
  if (nrow(new_annotation) == 0) {
    stop("No currently annotated site matches with your subset!")
  }
  new_annotation <- .generate_name_positions(new_annotation, "rna", "rnapos")
  ribo <- annotate_site(ribo,new_annotation)
  
  return(ribo)
  
}