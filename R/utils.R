#' Check if a vector of strings correspond to valid metadata in a given RiboClass
#'
#' This function returns nothing if each string matches with a metadata.
#' It will stop and display a formatted error message if one or more strings in
#' the vector are not valid metadata.
#' @param ribo A RiboClass.
#' @param metadata_name The vector of string to check against ribo's metadata.
#' @keywords internal
#'
check_metadata <- function(ribo,metadata_name) {
  # Get metadata columns that do not exist in ribo's metadata.
  unmatched_elts <- metadata_name[
    which(!(metadata_name %in% colnames(ribo[["metadata"]])))
    ]
  len_unmatched <- length(unmatched_elts)
  if(len_unmatched > 0) {
    
    cli::cli_abort(c(
      "You have supplied names that are not part of the RiboClass metadat.",
      "i" = "The RiboClass has the following metadata:
      {.val {colnames(ribo[['metadata']])}}.",
      "x" = "{len_unmatched} supplied name{?s} {?is/are} not part of the metadata:
      {.val {unmatched_elts}}."
    ))
    
  }
}
#' Check if a vector of strings correspond to valid samplenames in a given RiboClass
#'
#' This function returns nothing if each string matches with the samplenames
#' It will stop and display a formatted error message if one or more strings in
#' the vector are not valid samplenames
#' @param ribo A RiboClass.
#' @param sample_names The vector of string to check against ribo's samplenames
#' @keywords internal
#'
check_sample <- function(ribo, sample_names) {
  unmatched_elts <- sample_names[
    which(!(sample_names %in% names(ribo[["data"]])))
  ]
  len_unmatched <- length(unmatched_elts)
  
  if(len_unmatched > 0) {
    cli::cli_abort(c(
      "You have supplied samplenames that are not part of the supplied RiboClass",
      "i" = "Sample names in RiboClass : {.val {names(ribo[['data']])}}.",
      "x" = "{len_unmatched} supplied name{?s} {?is/are} not part of the RiboClass' samplenames:
      {.val {unmatched_elts}}."
    ))
  }
}