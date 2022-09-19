#' Create a RiboClass from count files and metadata.
#' 
#' @param count_path (required) path to the data folder containing count files
#' @param metadata  data frame or path to a CSV file containing metadata.
#' @param count_sep delimiter used in genomecov (for csv file only)
#' @param metadata_sep delimiter used in metadata (for csv file only)
#' @param count_header boolean, specify if count files have a header or not. 
#' @param count_value column containing count values
#' @param count_rnaid name or position of the column containing the name of the RNA in count data.
#' @param count_pos name or position of the column containing the site's position in count data.
#' @param metadata_key name or position of the column containing the samples' filename
#' @param metadata_id if it exists, name or position of the column containing the sample name
#'
#' @description 
#' Read ribomethseq count files and their associated metadata and turn them into a RiboClass.
#' __This is the entrypoint for the Riboscore package__, as all other functions use the RiboClass as an input.
#' 
#' 
#' @details
#' The RiboClass object is a S3 Class with three elements :
#' 
#' 
#'  __data__ : a list of dataframe, each corresponding to a sample.
#' Each dataframe starts with the following columns : 
#' * rna : the name of the RNA for a given position
#' * rnapos : the position __on the current RNA__
#' * count : the number of read starting (5'end) or ending (3'end) at this position
#' * site : the name of the site, which will be empty after the RiboClass creation. To fill this column, use \code{\link{annotate_site}}.
#' 
#' 
#' __metadata__ : a dataframe containing all information related to each sample.
#' This is particularly useful for technical biases correction (with \code{\link{adjust_bias}}) and plot functions (for coloring or add an extra layer of informations).
#' If no metadata are given, an empty one will be generated. You can also generate an empty one by yourself using \code{\link{generate_metadata_df}}.
#' 
#' __rna_names__ : a dataframe containing original and current rna names.
#' The name of RNAs can be renamed for the sake of clarity on plots (with \code{\link{rename_rna}}), but the original ones can contain important information, like NCBI ID.
#' 
#' 
#' the path given in __count_path__ should contains only necessary CSV files (one per sample).
#' While the directory structure is not important, make sure each sample has an unique filename.
#' 
#' The path to the csv file or the dataframe given in __metadata__ must contains a filename column, as this will serve to link metadata with the dataframes in data during the RiboClass creation.
#' 
#' 
#' @md
#' 
#' @example 
#' csv_files_path <- system.file("extdata", "miniglioma", package = "Riboscore")
#' 
#' @return a RiboClass
#' @export
create_riboclass <- function(count_path,
                             metadata = NULL,
                             count_sep = "\t",
                             metadata_sep = ",",
                             count_header = FALSE,
                             count_value = 3,
                             count_rnaid = 1,
                             count_pos = 2,
                             metadata_key = "filaname",
                             metadata_id = NULL) {
  
  # named_position => default_posname
  # col des positions connus => siteID
  
  # read count files
  #rna_counts_dt <- .read_count_files(count_path,count_sep,count_header,count_rnaid,count_pos,count_value)
  
  #create a table containing rna names
  
  #TODO rna_names -> rnaname
  #TODO has_cscore -> check premier tableau
  
  
  #loading metadata
  if(is.null(metadata)) {
    metadata <- generate_metadata_df(count_path,create_samplename_col = F)
    rna_counts_dt <- .read_count_files(count_path,count_sep,count_header,count_rnaid,count_pos,count_value)
    if(length(rna_counts_dt) == 0 ) stop("ERROR : no files were loaded")
    rna_names_df <- .generate_rna_names_table(rna_counts_dt[[1]])
  }
  
  else {
    
    if(is.character(metadata)) {
      metadata <- utils::read.csv(metadata, sep = metadata_sep)
      
    }
    
    if(!is.data.frame(metadata)) {
      stop("metadata must be a dataframe or a path to a csv file !")
    }
    
    names(metadata)[names(metadata) == metadata_id] <- "samplename"
    
    rownames(metadata) <- metadata[,"samplename"]
    
    # read count data
    rna_counts_dt <- .read_count_files(count_path,count_sep,count_header,
                                       count_rnaid,count_pos,count_value,
                                       files_to_keep = as.character(metadata[,metadata_key]))
    
    if(length(rna_counts_dt) == 0) stop(paste0("ERROR! No file has the filenames specified in your '",colnames(metadata[metadata_key]),"' column"))
    # generate RNA names table
    rna_names_df <- .generate_rna_names_table(rna_counts_dt[[1]])
    
    # Rename sample in counts list according to the names in metadata
    names(rna_counts_dt) <- metadata[,"samplename"][match(names(rna_counts_dt), metadata[,metadata_key])]
    
    # order samples by metadata
    
    rna_counts_dt <- rna_counts_dt[metadata[,"samplename"]]    
    
  }
  
  
  # Merge metadata and counts in the single named list.
  return_list <- list(data = rna_counts_dt, metadata = metadata,rna_names = rna_names_df, has_cscore = FALSE)
  class(return_list) <- "RiboClass"
  
  cat("[SUCCESS] Your data have been imported and the following RiboClass has been created :\n")
  print(return_list)
  
  return(return_list)
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
#' @param files_to_keep if specified, only the file following files_to_keep will be kept
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
           files_to_keep = NULL) {
    
    if(!dir.exists(path_to_files)) stop("the path given for the csv files does not exist or is not a directory !")
    
    rna_counts_fl <-
      list.files(path_to_files, recursive = T, full.names = T)
    
    # Check if there are files to keep
    # if yes than keep only the needed files
    if (!is.null(files_to_keep)) {
      pat <- paste0("\\b(", paste(files_to_keep, collapse="|"), ")\\b")
      rna_counts_fl <- rna_counts_fl[grep(pat, rna_counts_fl)]
    }
    
    # 1) Check if there is any duplicated name in the filename list
    
    if (anyDuplicated(basename(rna_counts_fl)) > 0) {
      stop("ERROR: some samples share the same filename!")
    }
    
    
    # 2) reorder cols for each count table and name them like this :
    # RNA | position_on_rna | count
    # and add a siteID column with a default value of NA
    rna_counts_dt <-
      lapply(rna_counts_fl, utils::read.csv, sep = sep, header = header)
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
      sapply(names(rna_counts_dt), function(x)
        check_sample_positions(
          sample_1 = rna_counts_dt[[1]],
          sample_1_name = reference_sample_name,
          sample_2 = rna_counts_dt[[x]],
          sample_2_name = x
        ))
    
    failing_samples <-
      names(sample_check_results[which(sample_check_results == F)])
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


#' Aggregate results into a single matrix
#'
#' For a given column in data, this function will generate a dataframe with all samples
#' and all positions (if only_annotated is false) or only positions with a siteID (if only_annotated is true).
#'
#' @param ribo a RiboClass object
#' @param col column in data you want extract data from (cscore or count).
#' @param position_to_rownames if true, position will be included as a rowname. They will in a new column otherwise.
#' @param only_annotated if true, return a dataframe with only annotated sites. Return all sites otherwise.
#' @return a dataframe 
#' @export
#'
#' @examples
#' data("ribo_toy")
#' count_df <- extract_data(ribo_toy,"count")
extract_data <- function(ribo, col = "cscore", position_to_rownames =F, only_annotated = F) {
  
  named_position <- NULL # NSE fix
  #TODO : sample_list -> ribo
  #TODO : check if position_to_rownames is useful ? -> trash
  #TODO : check if col exist
  #TODO : select only siteIdentified sites
  #The rows of this matrix correspond to the positions on the rRNA
  
  
   col <- tolower(col)
   
   if(!(col %in% colnames(ribo[["data"]][[1]]))) {
     
     stop(paste(col, "is not a column in the data !"))
   }
  
  
  sample_list <- ribo[["data"]]
  sample_list_nm <- names(sample_list)
  position_list <- sample_list[[3]][,"named_position"]
  matrix_all <- data.frame(named_position = position_list)
  
  for(sample_nm in sample_list_nm) {
    sample_df <- sample_list[[sample_nm]]
    sample_df <- sample_df[,c("named_position",col)]
    matrix_all <- dplyr::full_join(matrix_all,sample_df,by="named_position")
    names(matrix_all)[length(names(matrix_all))] <- sample_nm 
  }
  matrix_all <- matrix_all[match(position_list,matrix_all[,"named_position"]),]
  
  if(position_to_rownames){
    
    row.names(matrix_all) <- matrix_all[,"named_position"]
    matrix_all <- subset(matrix_all,select = -named_position)
  }
  
  return(matrix_all)
  
  
}

#' Rename RNAs in your RiboClass
#'
#'  
#' 
#' @param ribo a RiboClass object
#' @param new_names the new names for your RNA (by order of rna size)
#'
#' 
#'
#' @return a RiboClass with updated rna names in data.
#' @export
#'
#' @examples
#' ribo_toy <- rename_rna(ribo_toy ,c("5S","5.8S","18S","28S"))

rename_rna <- function(ribo,new_names=c("5S","5.8S","18S","28S")) {
  
  sample_list <- ribo[["data"]]
  rna_names <- ribo[["rna_names"]]
  rna_names[3] <- new_names
  
  #change the RNA names inside each sample
  sample_list_renamed <- lapply(sample_list, function(x) {
    x[,1] <- rna_names[,3][match(x[,1], rna_names[,2])]
    return(x)
  })
  rna_names[2] <- rna_names[3]
  rna_names <- rna_names[,1:2]
  # Update nomenclature according to the new RNAs
  sample_list_renamed <- .generate_riboclass_named_position(sample_list_renamed,1,2)

  ribo[["data"]] <- sample_list_renamed
  ribo[["rna_names"]] <- rna_names
  return(ribo)
}

#' (internal) Generate a table with former and current rna names for a given RiboClass
#' @param count_df a count file dataframe, containing a columns with rna names
#'
#' @return dataframe with former and current rna names
#' 
#'
#' @keywords internal
.generate_rna_names_table <- function(count_df) {
  
   rna_counts <- as.data.frame(sort(table(count_df[,1])))

  rna_names_df <- data.frame(original_name = rna_counts[[1]], current_name =rna_counts[[1]])
  return(rna_names_df)
}

#' Generate the default name for positions
#' This function is used to generate the named_position column in a given sample.
#' It is always called when generating a RiboClass.
#' @param sample_count_list list of count dataframe (one sample = one count dataframe)
#' @param rna_col name or position of the column containing RNA names
#' @param rnapos_col name or position of the column containing position in RNA
#' @keywords internal
#' @return a list of sample data with the added "named_position" column
#' 
#'
.generate_riboclass_named_position <- function(sample_count_list,rna_col,rnapos_col) {
  #combine the RNA name and the position on this RNA to form the row names.
  sample_count_list_named <- lapply(sample_count_list,function(x){ 
    x["named_position"] <-  paste(x[,rna_col], formatC(x[,rnapos_col],width = 4,flag = "0"), sep = "_"); x 
  })
  return(sample_count_list_named)
}


#' Update count values inside a RiboClass with a matrix of values
#' 
#' @param ribo a RiboClass object
#' @param update_matrix a position x sample matrix containing the new values 
#'
#' @return a RiboClass with updated values
#' @export
#'
update_ribo_count_with_matrix <- function(ribo, update_matrix) {
  #first, check if we have the sample name in our column
  update_df <- as.data.frame(update_matrix)
  col_names <- sort(names(update_df))
  riboclass_names <- sort(names(ribo[["data"]]))
  
  if(!identical(col_names,riboclass_names)) {
    stop("mismatch between samplenames and matrix's samples names")
  }
  
  # For each sample in the RiboClass, replace count values with matrix's ones
  
  count_list <- ribo[["data"]]
  
  for(sample in names(count_list)) {
    count_list[[sample]]["count"] <- update_df[sample]
  }
  
  ribo[["data"]] <- count_list
  
  return(ribo)
  
  
}




#' Regroup samples by condition and calculate mean for each condition
#' 
#' @description An helper function that will give the mean of the cscore or count of all samples by condition, for each position. The standard deviation is also given.
#' this can be used to create boxplot with ggplot.
#'
#' @param ribo a RiboClass object
#' @param metadata_condition name or position of the column __in metadata__ containing the condition
#' @param value name or position of the column containing the values on which mean by condition is calculated.
#' 
#' @importFrom dplyr %>%
#' @importFrom rlang sym
#' @return a dataframe with the mean for each condition for a selected value
#' @export
#' @md
#' @example 
#' mean_df <- mean_samples_by_conditon(ribo_toy,"count","condition")
#' 
mean_samples_by_conditon <- function(ribo,value, metadata_condition) {
  
  named_position <- NULL # NSE fix
  ribo_list <- ribo[["data"]]
  ribo_names <- names(ribo_list)
  ribo_list_named <- lapply(ribo_names, function(x){
    ribo_list[[x]]["sample"] <- x
    return(ribo_list[[x]])
  })
  
  ribo_concat <- dplyr::bind_rows(ribo_list_named)
  
  metadata <- ribo[["metadata"]]
  ribo_concat[metadata_condition] <- metadata[,metadata_condition][match(ribo_concat[,"sample"], metadata[,"samplename"])]
  ribo_condition <- ribo_concat %>% dplyr::group_by(named_position, !!sym(metadata_condition)) %>% dplyr::summarise(mean = mean(!!sym(value)), sd = stats::sd(!!sym(value)))
  return(ribo_condition)
}


#' Remove a RNA among all samples
#' 
#' @param ribo a RiboClass
#' @param name_rna_to_remove name of the rna to remove
#'
#' @return a RiboClass without the specified RNA
#' @export
#'
ribo_remove_rna <- function(ribo, name_rna_to_remove) {
  #TODO : use lapply instead
  for(name in names(ribo$counts)) {
    
    df <- ribo$counts[[name]]
    df <- df[df["RNA"] != name_rna_to_remove,]
    ribo$counts[[name]] <- df
    
  }
  
  
  
  return(ribo)
  
}

#' Display RNA names
#'
#' @param ribo a RiboClass object
#'
#' @return
#' A vector with actual RNA names
#' @export
#'
#' @examples
#' show_RNA_names(ribo = ribo_toy)
show_RNA_names <- function(ribo = NULL) {
  if(is.null(ribo)) {stop("A ribo class object should be provided")}
  RNA_names <- ribo[["rna_names"]][[2]]
  return(RNA_names)
}

#' print() to display basic informations about a RiboClass
#'
#' @param x a RiboClass
#' @param ... base print params
#'
#' @return
#' @export
#' @keywords internal
#' @examples
print.RiboClass <- function(x,...){ 
  cat(as.character(paste("a RiboClass with", length(x[["data"]]),"samples and", nrow(x[["rna_names"]]), "RNA(s) :")))
   rna <- table(as.factor(x[["data"]][[1]][["rna"]]))
   rna_df <- as.data.frame(rna)
   rna_df <- rna_df[order(rna_df[,2]),]
   cat(paste0("\nName : ",rna_df[,1] , ", length : ", rna_df[,2]))
}