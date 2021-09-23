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