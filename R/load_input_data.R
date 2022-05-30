#' Read ribomethseq count files and their associated metadata and turn them into a riboclass.
#' This should be considered as the entrypoint for the RMS package, as all other functions use the riboclass as an input.
#' 
#' @param counts_path: path to the data folder with the count files
#' @param metadata: data frame or path to the CSV file containing metadata
#' @param columns_names : list of names to set for each column (from left to right)
#' @param count_sep: delimiter used in genomecov csv
#' @param metadata_sep: delimiter used in metadata csv
#' @param count_header: boolean, specify if count files have a header or not. 
#' @param count_rna: name or position of the column containing the name of the RNA in counts data
#' @param count_rnapos: name or position of the column containing the position on an RNA in counts data
#' @param metadata_filename: name or position of the column containing the filename
#' @param metadata_samplename: name or position of the column containing sample name
#' @return a riboclass
#' @export
create_riboclass <- function(counts_path,
                             metadata = NA,
                             count_sep = "\t",
                             metadata_sep = ";",
                             count_header = FALSE,
                             count_col = 3,
                             count_rna = 1,
                             count_rnapos = 2,
                             metadata_filename = 1,
                             metadata_samplename = 2) {
  
  # named_position => default_posname
  # col des positions connus => siteID
  
  # read count files
  #rna_counts_dt <- read_count_files(counts_path,count_sep,count_header,count_rna,count_rnapos,count_col)
  
  #create a table containing rna names
  
  #TODO rna_names -> rnaname
  #TODO has_cscore -> check premier tableau

  
  #loading metadata
  if(!is.data.frame(metadata)) {
    
    if(is.character(metadata)) {
      metadata <- read.csv(metadata, sep = metadata_sep)
      rownames(metadata) <- metadata[,"samplename"]
    }
    else if (is.na(metadata)) {
      
      metadata <- generate_metadata_df(counts_path,create_samplename_col = F)
      rna_counts_dt <- read_count_files(counts_path,count_sep,count_header,count_rna,count_rnapos,count_col)
      rna_names_df <- generate_rna_names_table(rna_counts_dt[[1]])
      
    }
    

  }
  else if(is.data.frame(metadata)) { 
  
    names(metadata)[names(metadata) == metadata_samplename] <- "samplename"
    
    rownames(metadata) <- metadata[,"samplename"]
    
    # read count data
    rna_counts_dt <- read_count_files(counts_path,count_sep,count_header,
                                      count_rna,count_rnapos,count_col,
                                      files_to_keep = as.character(metadata[,metadata_filename]))
    # generate RNA names table
    rna_names_df <- generate_rna_names_table(rna_counts_dt[[1]])
    
    # Rename sample in counts list according to the names in metadata
    names(rna_counts_dt) <- metadata[,"samplename"][match(names(rna_counts_dt), metadata[,metadata_filename])]
    
    # order samples by metadata
    
    rna_counts_dt <- rna_counts_dt[metadata[,"samplename"]]
  }
  
  else {
    stop("incorrect metadata given")
  }
  
  # Merge metadata and counts in the single named list.
  return_list <- list(data = rna_counts_dt, metadata = metadata,rna_names = rna_names_df, has_cscore = FALSE)
  class(return_list) <- "RiboClass" # TODO : create a real constructor
  return(return_list)
}



#' Generate a metadata dataframe, given a list of files
#' @param counts_folder_path: the path where count files are stored
#' @param create_samplename_col: generate a sample name col with filename by default
#' @param stop_symbol: keep the filename until the stop symbol is reached
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
#'
#' @return A boolean indicating if the two samples are identical.
#' 
#'
#' @examples
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
