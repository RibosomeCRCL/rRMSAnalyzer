#' Load csv files from GenomeCov and their associated metadata. Create a RiboClass.
#' @inheritParams new_riboclass
#' @inheritParams compute_cscore
#' @return a RiboClass
#' @export
#' 
#' @description Import your count CSV files and the metadata to create a RiboClass.
#' The RiboClass is used by rRMSAnalyzer package for all analyses. 
#' 
#' __This function serves as the entrypoint of rRMSAnalyzer.__
#' 
#' @md
#' @details
#' load_ribodata is a wrapper of \code{\link{new_riboclass}} and \code{\link{compute_cscore}}.
#' @seealso new_riboclass  
#'
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
  
  ribo <- new_riboclass(count_path,
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
  cli::cli_alert_success("Your data has been successfully loaded!")
  return(ribo)
  
}