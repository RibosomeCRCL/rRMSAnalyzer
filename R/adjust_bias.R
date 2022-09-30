#' Correct batch effect with ComBatseq method
#' 
#' @param ribo a RiboClass object, see constructor : \code{\link{create_riboclass}}
#' @param batch name of the column in metadata that contains the batch number.
#' @param ... Parameters for sva's \code{\link[sva]{ComBat_seq}} function.
#' @return RiboClass with adjusted values and recomputed cscore.
#' @export
#' 
#' @details 
#' You must have a column with the batch number for each sample in your RiboClass’s metadata.
#'  
#' @references Yuqing Zhang, Giovanni Parmigiani, W Evan Johnson, ComBat-seq: batch effect adjustment for RNA-seq count data, NAR Genomics and Bioinformatics, Volume 2, Issue 3, 1 September 2020, lqaa078, https://doi.org/10.1093/nargab/lqaa078
#'
#'
#' @examples
#' data("ribo_toy")
#' ribo_toy_two <- keep_ribo_samples(ribo_toy,c("4283","RNA1","4272","RNA2"))
#' ribo_toy_adjusted <- adjust_bias(ribo_toy_two,"run") 
#' 
adjust_bias <- function(ribo, batch,...) {
  # Get the count matrix.
  # example :
  #            VMT11.csv VMT12.csv VMT13.csv VMT14.csv VMT15.csv VMT16.csv
  # 18S_1          1         1         1         1         1         1
  # 18S_2        221       344       285       486       444       504
  # 18S_3         84       119       123       172       132       190
  # 18S_4         44        68        65       114        99       116
  # 18S_5         84       103       154       262       218       240
  # 18S_6        100       137       194       336       247       264
  
  if(batch %in% colnames(ribo[["metadata"]])) {
    matrix_ribo <- extract_data(ribo,"count",position_to_rownames = T)
    #reorganize column according to metadata and convert DF to matrix (otherwise, ComBat_seq won't work)
    matrix_ribo <- as.matrix(matrix_ribo[,c(ribo[["metadata"]][["samplename"]])])
    adjusted_matrix <- sva::ComBat_seq(matrix_ribo,batch = ribo[["metadata"]][[batch]],...)
    
    ribo_updated <- .update_ribo_count_with_matrix(ribo,adjusted_matrix)
    
    if(ribo_updated[["has_cscore"]]) {
      cat("Recomputing c-score with the following parameters :\n")
      cat(paste("c-score method :", ribo_updated[["cscore_method"]],"\n"))
      cat(paste("flanking window :", ribo_updated[["cscore_window"]],"\n"))
      
      ribo_updated <- compute_cscore(ribo_updated, ribo_updated[["cscore_window"]],ribo_updated[["cscore_method"]],2)
    }
    

    
    ribo_updated["combatSeq_count"] <- TRUE
    ribo_updated["col_used_combatSeq"] <- batch
    return(ribo_updated)
  }
  else {
    stop(paste0("\"",batch, "\" is not a valid metadata."))
  }

  
}

