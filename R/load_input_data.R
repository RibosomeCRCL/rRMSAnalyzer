#' Load csv files from GenomeCov and their associated metadata. Create a RiboClass.
#' @inheritParams create_riboclass
#' @inheritParams compute_cscore
#' @return
#' @export
#' 
#' @description Riboscore's entrypoint. 
#'
#' @examples
load_ribodata <- function(count_path,
                          metadata = NULL,
                          count_sep = "\t",
                          metadata_sep = ",",
                          count_header = FALSE,
                          count_col = 3,
                          count_rna = 1,
                          count_rnapos = 2,
                          metadata_filename = 1,
                          metadata_samplename = 2,
                          flanking=6,
                          method = "median",
                          ncores = 1) {
  
  ribo <- create_riboclass(count_path,
                           metadata,
                           count_sep,
                           metadata_sep,
                           count_header,
                           count_col,
                           count_rna,
                           count_rnapos,
                           metadata_filename,
                           metadata_samplename)
  
  ribo <- compute_cscore(ribo,flanking,method,ncores)
  
  return(ribo)
  
}



#' Generate a metadata dataframe, given a path.
#' 
#' @description 
#' Generate an "empty" dataframe with the required filename and samplename columns to include in a riboclass.
#' You are free to add any other columns for your metadata.
#' If you do not use stop_symbol, filename and samplename will be identical. Feel free to modify the latter.
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
#' @export
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



#' Check if two samples share the same positions
#'
#' @param sample_1 First sample to check against sample_2
#' @param sample_2 Second sample to check against sample_1
#' @param sample_1_name name of sample 1
#' @param sample_2_name name of sample 2
#' @keywords internal
#' @return A boolean indicating if the two samples are identical.
check_sample_positions <- function(sample_1, sample_2,sample_1_name,sample_2_name) {
  sample_size <- length(sample_1[, "named_position"])

  if (length(sample_2[, "named_position"]) == sample_size) {
    if (sum(sample_1["named_position"] == sample_2["named_position"]) == sample_size)
      return(TRUE)
    else {
      errored_positions_file_1 <- setdiff(sample_1[["named_position"]],sample_2[["named_position"]])
      errored_positions_file_2 <- setdiff(sample_2[["named_position"]],sample_1[["named_position"]])
      paste("[WARNING] Two samples have differents positions! \n Mismatch at ", sample_1_name, ": ",paste(errored_positions_file_1,collapse = ", "),"\n Mismatch at ",sample_2_name, ": ",paste(errored_positions_file_2,collapse = ", "))
      return(FALSE)
    }
    
  }
  else {
    warning(paste("[WARNING]", sample_1_name, "and", sample_2_name ," have different size!"))
    return(FALSE)
  }
  
}
