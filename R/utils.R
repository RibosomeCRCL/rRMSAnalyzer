#' Generate "named_position" column
#'
#' @param df dataframe with both rna and pos columns
#' @param rna name or position of rna column
#' @param pos name or position of position in rna column
#' @keywords internal
#' @return a dataframe
generate_name_positions <- function(df, rna, pos) {
  df["named_position"] <- paste(df[,rna], formatC(df[,pos],width = 4,flag = "0"), sep = "_")
  return(df)
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