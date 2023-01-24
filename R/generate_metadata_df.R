#' Generate a metadata dataframe, given a path.
#' 
#' @description 
#' Generate an "empty" dataframe with the required filename and samplename columns to include in a riboclass.
#' You are free to add any other columns for your metadata.
#' If you do not use stop_symbol, filename and samplename will be identical.
#' Feel free to modify the latter.
#' 
#' __Do not modify the filename column, unless you have changed the filenames on disk.__ Otherwise, it would prevent \code{\link{create_riboclass}} from linking data and metadata...
#' 
#' @details 
#' stop_symbol is an "helper" parameter to automatically trim your filename into samplename.
#' 
#' @md 
#' @param counts_folder_path the path where count files are stored
#' @param create_samplename_col generate a sample name col with filename by default
#' @param stop_symbol keep the filename until the stop symbol is reached.
#' 
#' @keywords internal
generate_metadata_df <- function(counts_folder_path,
                                 create_samplename_col=T,
                                 stop_symbol=NA) {
  
  sample_filenames <- basename(list.files(counts_folder_path, recursive = T))
  
  if(create_samplename_col) {
    sample_name <- sample_filenames
    #if a symbol has been given to shorten the name
    if(!is.na(stop_symbol)) {
      sample_name <- stringr::str_extract(sample_name,
                                          paste0("^([^", stop_symbol,"])+"))
    }
    metadata_template <- data.frame(filename = sample_filenames, samplename=sample_name)
  }
  else {
    metadata_template <- data.frame(filename = sample_filenames)
  }
  
  return(metadata_template)
}
