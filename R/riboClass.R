#' Import and transform a list of count files for a riboClass
#' @description This function should not be used by itself. It is better to generate a full riboClass.
#' @param path_to_files path to the folder containing count files
#' @param sep CSV seperator 
#' @param rna_col position of the column containing RNA name
#' @param position_col position of the column containing the genomic position
#' @param count_col position of the column contaning the count for the given genomic position 
#'
#' @return
#'
#' @examples
read_count_files <-
  function(path_to_files,
           sep,
           header,
           rna_col,
           position_col,
           count_col) {
    rna_counts_fl <-
      list.files(path_to_files, recursive = T, full.names = T)
    
    # 1) Check if there is any duplicated name in the filename list
    
    if (anyDuplicated(basename(rna_counts_fl)) > 0) {
      stop("ERROR: some samples share the same filename!")
    }
    
    
    # 2) reorder cols for each count table and name them like this :
    # RNA | position_on_rna | count
    rna_counts_dt <-
      lapply(rna_counts_fl, read.csv, sep = sep, header = header)
    rna_counts_dt <- lapply(rna_counts_dt, function(x) {
      x <- x[, c(rna_col, position_col, count_col)]
      colnames(x) <- c("RNA", "Position_on_RNA", "Count")
      return(x)
    })
    
    # 3) combine both the RNA and genomic position columns to generate the column "named_position".
    rna_counts_dt <-
      generate_riboclass_named_position(rna_counts_dt, "RNA", "Position_on_RNA")
    
   
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
      print("[SUCCESS] all samples have identical positions!")
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
#' @param sample_list 
#' @param col_to_keep 
#' @param position_to_rownames 
#'
#' @return
#' @export
#'
#' @examples
aggregate_samples_by_col <- function(sample_list, col_to_keep, position_to_rownames =F) {
  
  #The rows of this matrix correspond to the positions on the rRNA
  sample_list_nm <- names(sample_list)
  position_list <- sample_list[[3]][,"named_position"]
  matrix_all <- data.frame(named_position = position_list)
  
  for(sample_nm in sample_list_nm) {
    sample_df <- sample_list[[sample_nm]]
    sample_df <- sample_df[,c("named_position",col_to_keep)]
    matrix_all <- full_join(matrix_all,sample_df,by="named_position")
    names(matrix_all)[length(names(matrix_all))] <- sample_nm 
  }
  matrix_all <- matrix_all[match(position_list,matrix_all[,"named_position"]),]
  
  if(position_to_rownames){
    
    row.names(matrix_all) <- matrix_all[,"named_position"]
    matrix_all <- subset(matrix_all,select = -named_position)
  }
  
  return(matrix_all)
  
  
}

#' Update rna names
#'
#' @param sample_list
#' @param rna_col 
#' @param rna_names 
#'
#' @return
#' @export
#'
#' @examples
#' 
#TODO : add "actual rna names" col
update_riboclass_rna_names <- function(ribo) {
  
  sample_list <- ribo[["counts"]]
  rna_names <- ribo[["rna_names"]]
  
  sample_list_renamed <- lapply(sample_list, function(x) {
    x[,1] <- rna_names[,3][match(x[,1], rna_names[,2])]
    return(x)
  })
  
  sample_list_renamed <- generate_riboclass_named_position(sample_list_renamed,1,2)
  rna_names[,2] <- rna_names[,3]
  rna_names[,3] <- NA
  ribo[["counts"]] <- sample_list_renamed
  ribo[["rna_names"]] <- rna_names
  return(ribo)
}

#' Generate a table with current rna names for a given riboClass
#' The table also contains 
#' @param count_df 
#'
#' @return
#' @export
#'
#' @examples
generate_rna_names_table <- function(count_df) {
  rna_names <- unique(count_df[,1])
  rna_names_df <- data.frame(original_name = rna_names, current_name =rna_names,new_rna_name = NA )
  return(rna_names_df)
}

#' Generate the default name for positions
#' This function is used to generate the named_position column in a given sample.
#' It is always called when generating a riboClass.
#' @param sample_count_list 
#' @param rna_col 
#' @param rnapos_col 
#'
#' @return
#' @export
#'
#' @examples
generate_riboclass_named_position <- function(sample_count_list,rna_col,rnapos_col) {
  #combine the RNA name and the position on this RNA to form the row names.
  sample_count_list_named <- lapply(sample_count_list,function(x){ 
    x["named_position"] <-  paste(x[,rna_col], x[,rnapos_col], sep = "_"); x 
  })
  return(sample_count_list_named)
}


#' Update count values inside a riboclass with a matrix of values
#' 
#' @param ribo a riboclass object
#' @param update_matrix a position x sample matrix containing the new values 
#'
#' @return a riboclass with updated values
#' @export
#'
#' @examples
update_ribo_count_with_matrix <- function(ribo, update_matrix) {
  
  #first, check if we have the sample name in our column
  update_df <- as.data.frame(update_matrix)
  col_names <- sort(names(update_df))
  riboclass_names <- sort(names(ribo[["counts"]]))
  
  if(col_names != riboclass_names) {
    error("mismatch between samplenames and matrix's samples names")
  }
  
  # For each sample in the riboclass, replace count values with matrix's ones
  
  count_list <- ribo[["counts"]]
  
  for(sample in names(count_list)) {
    count_list[[sample]]["Count"] <- update_df[sample]
  }
  
  ribo[["counts"]] <- count_list
  
  return(ribo)
  
  
}




#' Regroup samples by condition and calculate mean for each condition
#'
#' @param ribo a riboClass object
#' @param metadata_condition condition to group samples by
#' @param value value to calculate mean by condition
#'
#' @return
#' @export
#'
#' @examples
mean_samples_by_conditon <- function(ribo,value, metadata_condition) {
  
  ribo_list <- ribo[["counts"]]
  ribo_names <- names(ribo_list)
  ribo_list_named <- lapply(ribo_names, function(x){
    ribo_list[[x]]["sample"] <- x
    return(ribo_list[[x]])
  })
  
  ribo_concat <- dplyr::bind_rows(ribo_list_named)
  
  metadata <- ribo[["metadata"]]
  ribo_concat[metadata_condition] <- metadata[,metadata_condition][match(ribo_concat[,"sample"], metadata[,"samplename"])]
  ribo_condition <- ribo_concat %>% group_by(named_position, !!sym(metadata_condition)) %>% summarise(mean = mean(!!sym(value)), sd = sd(!!sym(value)))
  return(ribo_condition)
}


#' Display RNA names
#'
#' @param ribo a riboClass object
#'
#' @return
#' A vector with actual RNA names
#' @export
#'
#' @examples
#' show_RNA_names(ribo = ribo)
show_RNA_names <- function(ribo = NULL) {
  if(is.null(ribo)) {stop("A ribo class object should be provided")}
  RNA_names <- ribo[["rna_names"]][[2]]
  return(RNA_names)
}


