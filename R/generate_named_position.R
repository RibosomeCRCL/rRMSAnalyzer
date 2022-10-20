#' Generate the default name for positions
#' This function is used to generate the named_position column in a given sample.
#' It is always called when generating a RiboClass.
#' @param sample_count_list list of count dataframe (one sample = one count dataframe)
#' @param rna_col name or index of the column containing RNA names
#' @param rnapos_col name or index of the column containing position in RNA
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

# (internal) generate named_position column using the 3 columns in the data.
generate_name_positions <- function(df, rna, pos) {
  df["named_position"] <- paste(df[,rna], formatC(df[,pos],width = 4,flag = "0"), sep = "_")
  return(df)
}