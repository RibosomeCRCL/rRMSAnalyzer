#' Load csv files from GenomeCov and their associated metadata. Create a RiboClass.
#' @inheritParams create_riboclass
#' @inheritParams compute_cscore
#' @return a RiboClass
#' @export
#' 
#' @description Import your count CSV files and the metadata to create a RiboClass.
#' The RiboClass is used by Riboscore package for all analyses. 
#' 
#' @details
#' load_ribodata is a wrapper of \code{\link{create_riboclass}} and \code{\link{compute_cscore}}.
#'
#' @examples
load_ribodata <- function(count_path,
                          metadata = NULL,
                          count_sep = "\t",
                          metadata_sep = ",",
                          count_header = FALSE,
                          count_value = 3,
                          count_rnaid = 1,
                          count_pos = 2,
                          metadata_key = "filename",
                          metadata_id = NULL,
                          flanking=6,
                          method = "median",
                          ncores = 1) {
  
  ribo <- create_riboclass(count_path,
                           metadata,
                           count_sep,
                           metadata_sep,
                           count_header,
                           count_value,
                           count_rnaid,
                           count_pos,
                           metadata_key,
                           metadata_id)
  
  ribo <- compute_cscore(ribo,flanking,method,ncores)
  
  return(ribo)
  
}