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
  
  sample_list <- ribo[["raw_counts"]]
  rna_names <- ribo[["rna_names"]]
  
  sample_list_renamed <- lapply(sample_list, function(x) {
    x[,1] <- rna_names[,2][match(x[,1], rna_names[,1])]
    return(x)
  })
  
  sample_list_renamed <- generate_riboclass_named_position(sample_list_renamed,1,2)
  
  ribo[["raw_counts"]] <- sample_list_renamed
  return(ribo)
}

#' Title
#'
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
  
