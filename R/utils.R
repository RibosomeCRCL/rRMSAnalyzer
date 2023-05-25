#' Check if a vector of strings correspond to valid metadata in a given RiboClass
#'
#' This function returns nothing if each string matches with a metadata.
#' It will stop and display a formatted error message if one or more strings in
#' the vector are not valid metadata.
#' @param ribo A RiboClass.
#' @param metadata_name The vector of string to check against ribo's metadata.
#'
check_metadata <- function(ribo,metadata_name) {
  # Get metadata columns that do not exist in ribo's metadata.
  unmatched_elts <- metadata_name[
    which(!(metadata_name %in% colnames(ribo[["metadata"]])))
    ]
  len_unmatched <- length(unmatched_elts)
  if(len_unmatched > 0) {
    
    cli::cli_abort(c(
      "You have supplied names that are not part of the RiboClass metadata",
      "i" = "The RiboClass has the following metadata:
      {.val {colnames(ribo[['metadata']])}}",
      "x" = "{len_unmatched} supplied name{?s} {?is/are} not part of the metadata:
      {.val {unmatched_elts}}"
    ))
    
  }
}