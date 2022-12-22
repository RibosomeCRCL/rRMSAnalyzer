# (internal) Generate named_position column using the 3 columns in the data.:
#  rna : index of the column containing the name of the RNA
#  pos : index of the column containing the (numbered) position
# width specifies the number width. With a width of 4, "1" is written as "OOO1".
#
# It returns a list of named position in the following format : rna_position
# Example : 18S_0023
.generate_name_positions <- function(df, rna, pos, width = 4) {
  df["named_position"] <- paste(df[,rna], formatC(df[,pos],width = width,flag = "0"), sep = "_")
  return(df)
}

# Generate the default name for positions
# This function is used to generate the named_position column in a given sample.
#  sample_count_list list of count dataframe (one sample = one count dataframe)
#  rna_col name or index of the column containing RNA names
#  rnapos_col name or index of the column containing position in RNA
# 
# Returns a list of sample data with the added "named_position" column
.generate_riboclass_named_position <- function(sample_count_list,rna_col,rnapos_col) {
  #combine the RNA name and the position on this RNA to form the row names.
  sample_count_list_named <- lapply(sample_count_list,function(x){ 
    x <- .generate_name_positions(x, rna_col,rnapos_col);x
  })
  return(sample_count_list_named)
}

