#' Create a RiboClass from count files and metadata.
#' 
#' @param count_path (required) path to the data folder containing count files.
#' @param metadata  Data frame or path to a CSV file containing metadata.
#' @param count_sep Delimiter used for the CSV files.
#' @param metadata_sep Delimiter used in metadata (if imported from file).
#' @param count_header Boolean, specify if count files have a header or not.
#' @param count_value Name or index of the column containing count values.
#' @param count_rnaid Name or index of the column containing the name of the RNA in count data.
#' @param count_pos Name or index of the column containing the site's position in count data.
#' @param metadata_key (required when metadata not null) Name or index of the column containing the samples' filename.
#' @param metadata_id Name or index of the column containing the sample name.
#'
#' @description 
#' Read ribomethseq count files and their associated metadata and merge them into a RiboClass.
#' __This is the constructor of the RiboClass S3 class.__ All other functions use the RiboClass as an input.
#' 
#' The constructor cannot be called directly. To create a RiboClass, please use \code{\link{load_ribodata}}.
#' 
#' @details
#' The RiboClass object is a S3 Class with three elements :
#' 
#' 
#'  __data__ : a list of dataframe, each corresponding to a sample.
#' Each dataframe starts with the following columns : 
#' * rna : the name of the RNA for a given position.
#' * rnapos : the position __on the current RNA__.
#' * count : the number of read starting (5'end) or ending (3'end) at this position.
#' * site : the name of the site, which will be empty after the RiboClass creation. To fill this column, use \code{\link{annotate_site}}.
#' 
#' 
#' __metadata__ : a dataframe containing all information related to each sample.
#' This is particularly useful for technical biases correction (with \code{\link{adjust_bias}}) and plot functions (for coloring or add an extra layer of informations).
#' If no metadata are given, an empty one will be generated.
#' 
#' __rna_names__ : a dataframe containing original and current rna names.
#' The name of RNAs can be renamed for the sake of clarity on plots (with \code{\link{rename_rna}}), but the original ones can contain important information, like NCBI ID.
#' 
#' 
#' the path given in __count_path__ should contains only necessary CSV files (one per sample).
#' While the directory structure is not important, make sure each sample has an unique filename.
#' 
#' The path to the csv file or the dataframe given in __metadata__ must contains a filename column, as this will serve to link metadata with the dataframes in data during the RiboClass creation.
#' @md
#' 
#' @return A RiboClass.
#'  
new_riboclass <- function(count_path,
                             metadata = NULL,
                             count_sep = "\t",
                             metadata_sep = ",",
                             count_header = FALSE,
                             count_value = 3,
                             count_rnaid = 1,
                             count_pos = 2,
                             metadata_key = "filename",
                             metadata_id = NULL) {
  
  
  #loading metadata
  if(is.null(metadata)) {
    metadata <- generate_metadata_df(count_path,create_samplename_col = FALSE)
    rna_counts_dt <- .read_count_files(count_path,count_sep,count_header,count_rnaid,count_pos,count_value)
    if(length(rna_counts_dt) == 0) stop("ERROR : no file was loaded")
    rna_names_df <- .generate_rna_names_table(rna_counts_dt[[1]])
  }
  
  else {
    
    if(is.character(metadata)) {
      if (file.exists(metadata)) {
        metadata <- utils::read.csv(metadata, sep = metadata_sep)
      }
      else {
        stop("the path specified for metadata does not exist or is not a file !")
      }
      
    }
    
    if(!is.data.frame(metadata)) {
      stop("metadata must be a dataframe or a path to a csv file !")
    }

    if(is.character(metadata_key) && !(metadata_key %in% names(metadata))) {
      stop(metadata_key, " (metadata_key param) is not a column in metadata")
    }

    if(is.character(metadata_id) && !(metadata_id %in% names(metadata))) {
      stop(metadata_id, " (metadata_id param) is not a column in metadata")
    }

    #rename the column specified in "metadata_id" to "samplename"
    if(is.character(metadata_id)) {
      names(metadata)[names(metadata) == metadata_id] <- "samplename"
    } else {
      names(metadata)[metadata_id] <- "samplename"
    }
    
    # Check if metadata has not duplicated filename or samplename
    if(anyDuplicated(metadata[metadata_key])) 
      stop("ERROR! Duplicated filename in metadata. Each sample must have an unique filename.")
    
    if(anyDuplicated(metadata["samplename"])) 
      stop("ERROR! Duplicated samplename in metadata. Each sample must have an unique samplename.")
    
    rownames(metadata) <- metadata[,"samplename"]
    
    # read count data
    rna_counts_dt <- .read_count_files(count_path,count_sep,count_header,
                                       count_rnaid,count_pos,count_value,
                                       metadata_filenames = as.character(metadata[,metadata_key]))
    
    
    # generate RNA names table
    rna_names_df <- .generate_rna_names_table(rna_counts_dt[[1]])
    
    # keep only metadatas that have associated data
    
    metadata <- metadata[which(metadata[,metadata_key] %in% names(rna_counts_dt)),]
    
    # Rename sample in counts list according to the names in metadata
    names(rna_counts_dt) <- metadata[,"samplename"][match(names(rna_counts_dt), metadata[,metadata_key])]
    
    # order samples by metadata and remove null samples
    rna_counts_dt <- rna_counts_dt[metadata[,"samplename"]]  
    rna_counts_dt <- rna_counts_dt[lengths(rna_counts_dt) != 0]
    
    
  }
  
  
  # Merge metadata and counts in RiboClass.
  RiboClass <- list(data = rna_counts_dt, metadata = metadata,rna_names = rna_names_df, has_cscore = FALSE)
  class(RiboClass) <- "RiboClass"
  
  return(RiboClass)
}


#' Import and transform a list of count files for a RiboClass
#' @description 
#' This internal function is used to generate the "data" part of the RiboClass.
#' It will read all CSV files and turn them into a list of dataframes.
#'
#' @param path_to_files path to the folder containing count files
#' @param sep CSV seperator 
#' @param rna_col position of the column containing RNA name
#' @param position_col position of the column containing the genomic position
#' @param count_value position of the column contaning the count for the given genomic position 
#' @param header does the files contain an header ?
#' @param metadata_filenames Filenames in metadata to check against
#'
#' @return a list of sample dataframes
#'
#' @keywords internal
.read_count_files <-
  function(path_to_files,
           sep,
           header,
           rna_col,
           position_col,
           count_value,
           metadata_filenames = NULL) {
    
    if(!dir.exists(path_to_files)) stop("the path given for the csv files does not exist or is not a directory !")
    
    rna_counts_fl <-
      list.files(path_to_files, recursive = TRUE, full.names = TRUE)
    
    # Check if the filenames on disk match filenames in metadata. Fail otherwise.
    if (!is.null(metadata_filenames)) {
      pat <- paste0("\\b(", paste(metadata_filenames, collapse="|"), ")\\b")
      rna_count_fl_old <- rna_counts_fl
      rna_counts_fl <- rna_counts_fl[grep(pat, rna_counts_fl)]
      
      if(length(metadata_filenames)> length(rna_counts_fl)) {
        rna_count_missing <- metadata_filenames[which(
          !is.element(metadata_filenames,basename(rna_counts_fl)))]
         warning(paste("File",rna_count_missing,"does not exist. Typo in metadata ?",collapse = "\n"))
      }
      
    }
    
    # 1) Check if there is any duplicated name in the filename list
    
    if (anyDuplicated(basename(rna_counts_fl)) > 0) {
      stop("ERROR: some samples share the same filename!")
    }
    
    # Check if there is a mismatch between the filenames in metadata 
    # and the filenames in the folder
    
    rna_counts_dt <-
      lapply(rna_counts_fl, utils::read.csv, sep = sep, header = header)
    
    # stop if total mismatch between filenames and the filename column given in metadata
    if(length(rna_counts_dt) == 0) stop(paste0("ERROR! No file has the filenames specified in the metadata column specified by metadata_key"))
    
    # check if there are less than 3 columns, which can happen when one fails to
    # specify the correct separator
    if(ncol(rna_counts_dt[[1]]) < 3) {
      stop("not enough columns in your count data !\n Check if you have specified the correct columns separator in count_sep")
    }
    
    # 2) reorder cols for each count table and name them like this :
    # RNA | position_on_rna | count
    # and add a siteID column with a default value of NA
    rna_counts_dt <- lapply(rna_counts_dt, function(x) {
      x <- x[, c(rna_col, position_col, count_value)]
      colnames(x) <- c("rna", "rnapos", "count")
      x["site"] <- NA
      return(x)
    })
    
    # 3) combine both the RNA and genomic position columns to generate the column "named_position".
    rna_counts_dt <-
      .generate_riboclass_named_position(rna_counts_dt, "rna", "rnapos")
    
    
    # 4) give a name for each element of the list
    #TODO: Check if the elements in RNA_counts_dt are in the same order as in RNA_counts_fl
    names(rna_counts_dt) <- basename(rna_counts_fl)
    
    # 5) using named_position, we check if samples share the same positions
    reference_sample_name <- names(rna_counts_dt)[1]
    
    sample_check_results <-
      vapply(names(rna_counts_dt), function(x) {
        .check_sample_positions(
          sample_1 = rna_counts_dt[[1]],
          sample_1_name = reference_sample_name,
          sample_2 = rna_counts_dt[[x]],
          sample_2_name = x
        )}, logical(1))
    
    failing_samples <-
      names(sample_check_results[which(sample_check_results == FALSE)])
    if (identical(failing_samples, character(0))) {
    }
    else {
      stop(paste(
        "[ERROR] The following samples have failed the positions check : ",
        paste(failing_samples, collapse = "; ")
      ))
    }
    
    
    
    return(rna_counts_dt)
  }
