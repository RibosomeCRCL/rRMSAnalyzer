#' Read ribomethseq count files and their associated metadata
#' 
#' @param counts_path: path to the data folder with the count files
#' @param metadata_path: path to the CSV file containing metadata
#' @param columns_names : list of names to set for each column (from left to right)
#' @param count_sep: delimiter used in genomecov csv
#' @param metadata_sep: delimiter used in metadata csv
#' @param counts_rna_col: name or position of the column containing the name of the RNA in counts data
#' @param counts_rnapos_col: name or position of the column containing the position on an RNA in counts data
#' @param metadata_filename_col: name or position of the column containing the filename
#' @param metadata_sample_name_col: name or position of the column containing sample name
#' @return a riboclass
#' @export
read_counts <- function(counts_path,
                        metadata_path=NA,
                        columns_names=c("RNA","Position_on_RNA","Count"),
                        count_sep="\t",
                        metadata_sep=";",
                        counts_rna_col = 1,
                        counts_rnapos_col = 2,
                        metadata_filename_col= 1,
                        metadata_sample_name_col = 2) {
  
  rna_counts_fl <- list.files(counts_path, recursive = T, full.names = T)

  # Check if there is any duplicated name in the list

  if(anyDuplicated(basename(rna_counts_fl)) > 0) {
    stop("ERROR: some samples share the same filename!")
  }
  
  rna_counts_dt <- lapply(rna_counts_fl, read.csv, sep = count_sep, header = F)
  rna_counts_dt <- lapply(rna_counts_dt, function(x){
   colnames(x) <- columns_names
   return(x)
  })
  #combine the RNA name and the position on this RNA to form the row names.
  rna_counts_dt <- lapply(rna_counts_dt,function(x){ 
    x["named_position"] <-  paste(x[,counts_rna_col], x[,counts_rnapos_col], sep = "_"); x 
    })
  #TODO: Check if the elements in RNA_counts_dt are in the same order as in RNA_counts_fl
  names(rna_counts_dt) = basename(rna_counts_fl)
  
  #Compare if each sample has the same samples as the reference
  #The first sample will be used as a reference.
  reference_sample_name <- names(rna_counts_dt)[1]
  
  # checked_samples <- lapply(names(rna_counts_dt),function(x) {
  #   check_sample_positions(rna_counts_dt[[x]],rna_counts_dt[[1]],x,reference_sample_name)
  # })
  
  sample_check_results <- sapply(names(rna_counts_dt),function(x) check_sample_positions(
                           sample_1 = rna_counts_dt[[1]],
                           sample_1_name= reference_sample_name,
                           sample_2 = rna_counts_dt[[x]],
                           sample_2_name = x))

  failing_samples <- names(sample_check_results[which(sample_check_results== F)])
  if(identical(failing_samples, character(0))) {
    print("[SUCCESS] all samples have identical positions!")
  }
  else {
    warning(paste("[WARNING] The following samples have failed the positions check : ", paste(failing_samples,collapse = "; ")))
  }
  
  
  #loading metadata
  if(is.na(metadata_path)) {
    metadata_df <- generate_metadata_df(counts_path,create_samplename_col = F)
  }
  else {
    metadata_df <- read.csv(metadata_path, sep = metadata_sep)
    # Rename sample in raw_counts according to the names in metadata
    names(rna_counts_dt) <- metadata_df[,metadata_sample_name_col][which(names(rna_counts_dt) == metadata_df[,metadata_filename_col])]
    
  }
  
  # Merge metadata and counts in the single named list.
  return_list <- list(raw_counts = rna_counts_dt, metadata = metadata_df)
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
    
    metadata_template <- data.frame(filename = sample_filenames, name=sample_name)
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
      errored_positions <- setdiff(sample_1[["named_position"]],sample_2[["named_position"]])
      errored_positions2 <- setdiff(sample_2[["named_position"]],sample_1[["named_position"]])
      warning(paste("[WARNING] Two samples have differents positions! \n Mismatch at ", sample_1_name, ": ",paste(errored_positions,collapse = ", "),"\n Mismatch at ",sample_2_name, ": ",paste(errored_positions2,collapse = ", ")))
      return(FALSE)
    }
    
  }
  else {
    warning(paste("[WARNING]", sample_1_name, "and", sample_2_name ," have different size!"))
    return(FALSE)
  }
  
}
