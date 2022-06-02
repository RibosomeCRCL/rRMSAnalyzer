
#' Correct batch effect with ComBat seq method
#'
#' @param ribo a riboclass
#' @param batch name of the column having batch effect for each sample (inside riboClass's metadata)
#'
#' @return a riboClass with adjusted counts
#' @export
#'
#' @examples 
adjust_bias <- function(ribo, batch) {
  # Get the count matrix.
  # example :
  #            VMT11.csv VMT12.csv VMT13.csv VMT14.csv VMT15.csv VMT16.csv
  # 18S_1          1         1         1         1         1         1
  # 18S_2        221       344       285       486       444       504
  # 18S_3         84       119       123       172       132       190
  # 18S_4         44        68        65       114        99       116
  # 18S_5         84       103       154       262       218       240
  # 18S_6        100       137       194       336       247       264
  
  matrix_ribo <- extract_data(ribo,"count",position_to_rownames = T)
  #reorganize column according to metadata and convert DF to matrix (otherwise, ComBat_seq won't work)
  matrix_ribo <- as.matrix(matrix_ribo[,c(ribo[["metadata"]][["samplename"]])])
  adjusted_matrix <- sva::ComBat_seq(matrix_ribo,batch = ribo[["metadata"]][[batch]])
  
  ribo_updated <- update_ribo_count_with_matrix(ribo,adjusted_matrix)
  ribo_updated["combatSeq_count"] <- TRUE
  ribo_updated["col_used_combatSeq"] <- batch
  return(ribo_updated)
  
}

