#' Create a riboClass from count files and metadata.
#' 
#' @param count_path path to the data folder containing count files
#' @param metadata data frame or path to a CSV file containing metadata
#' @param count_sep delimiter used in genomecov (csv file only)
#' @param metadata_sep delimiter used in metadata (csv file only)
#' @param count_header boolean, specify if count files have a header or not. 
#' @param count_col column containing count values
#' @param count_rna name or position of the column containing the name of the RNA in counts data.
#' @param count_rnapos name or position of the column containing the position on an RNA in counts data.
#' @param metadata_filename name or position of the column containing the filename
#' @param metadata_samplename name or position of the column containing the sample name
#'
#' @description 
#' Read ribomethseq count files and their associated metadata and turn them into a riboclass.
#' This should be considered as the entrypoint for the RMS package, as all other functions use the riboclass as an input.
#' 
#' 
#' 
#' @details 
#' Count data is expected to have the following the format (columns ordering is not important)
#' 
#'| RNA | Position_on_RNA | count |
#'|-----|-----------------|-------|
#'| 18S | 1245            | 30492  |
#'| 18S | 1246            | 19674  |
#'| 18S | 1247            | 23673  |
#' 
#' @md
#' 
#' If given, the Metadata table must have a **filename **and a **samplename** columns.
#' 
#' The riboClass will contain the follwing elements : 
#' 
#' 
#' @return a riboclass
#' @export
create_riboclass <- function(count_path,
                             metadata = NULL,
                             count_sep = "\t",
                             metadata_sep = ",",
                             count_header = FALSE,
                             count_col = 3,
                             count_rna = 1,
                             count_rnapos = 2,
                             metadata_filename = 1,
                             metadata_samplename = 2) {
  
  # named_position => default_posname
  # col des positions connus => siteID
  
  # read count files
  #rna_counts_dt <- .read_count_files(count_path,count_sep,count_header,count_rna,count_rnapos,count_col)
  
  #create a table containing rna names
  
  #TODO rna_names -> rnaname
  #TODO has_cscore -> check premier tableau

  
  #loading metadata
  if(is.null(metadata)) {
    metadata <- generate_metadata_df(count_path,create_samplename_col = F)
    rna_counts_dt <- .read_count_files(count_path,count_sep,count_header,count_rna,count_rnapos,count_col)
    rna_names_df <- .generate_rna_names_table(rna_counts_dt[[1]])
  }
  
  else {
    
    if(is.character(metadata)) {
      metadata <- utils::read.csv(metadata, sep = metadata_sep)

    }
    
    if(!is.data.frame(metadata)) {
      stop("metadata must be a dataframe or a path to a csv file !")
    }
    
    names(metadata)[names(metadata) == metadata_samplename] <- "samplename"
    
    rownames(metadata) <- metadata[,"samplename"]
    
    # read count data
    rna_counts_dt <- .read_count_files(count_path,count_sep,count_header,
                                      count_rna,count_rnapos,count_col,
                                      files_to_keep = as.character(metadata[,metadata_filename]))
    
    if(length(rna_counts_dt) == 0) stop(paste0("ERROR! No file has the filenames specified in your '",colnames(metadata[metadata_filename]),"' column"))
    # generate RNA names table
    rna_names_df <- .generate_rna_names_table(rna_counts_dt[[1]])
    
    # Rename sample in counts list according to the names in metadata
    names(rna_counts_dt) <- metadata[,"samplename"][match(names(rna_counts_dt), metadata[,metadata_filename])]
    
    # order samples by metadata
    
    rna_counts_dt <- rna_counts_dt[metadata[,"samplename"]]    
    
  }

  
  # Merge metadata and counts in the single named list.
  return_list <- list(data = rna_counts_dt, metadata = metadata,rna_names = rna_names_df, has_cscore = FALSE)
  class(return_list) <- "RiboClass" # TODO : create a real constructor
  return(return_list)
}



#' Generate a metadata dataframe, given a list of files
#' @param counts_folder_path the path where count files are stored
#' @param create_samplename_col generate a sample name col with filename by default
#' @param stop_symbol keep the filename until the stop symbol is reached
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
