#' Correct batch effect with ComBat-seq method
#' 
#' @description
#' Batch effect of RiboMethSeq data can be adjusted using the ComBat-seq method. adjust_bias is a wrapper to perform ComBat-seq adjustment. 
#' 
#' It will return a new RiboClass with adjusted read end count values and C-scores automatically recomputed with the same setup parameters.
#'
#' 
#' @param ribo a RiboClass object.
#' @param batch Name of the column in metadata that contains the batch number.
#' @param ncores Number of cores to use in case of multithreading.
#' @param ... Parameters to pass to sva's \code{\link[sva]{ComBat_seq}} function.
#' @return RiboClass with adjusted read end count values and automatically recomputed C-scores.
#' @export
#' 
#' @details 
#' You must have a column with the batch number for each sample in your RiboClass’s metadata.
#'  
#' @references Yuqing Zhang, Giovanni Parmigiani, W Evan Johnson, ComBat-seq: batch effect adjustment for RNA-seq count data, NAR Genomics and Bioinformatics, Volume 2, Issue 3, 1 September 2020, lqaa078, https://doi.org/10.1093/nargab/lqaa078
#'
#'
#' @examples
#' data('ribo_toy')
#' ribo_toy_two <- keep_ribo_samples(ribo_toy,c('S1','RNA1','S7','RNA2'))
#' ribo_toy_adjusted <- adjust_bias(ribo_toy_two,'run') 
 
adjust_bias <- function(ribo, batch,ncores = 1, ...) { # ... allows additional arguments to be passed to sva::ComBat_seq(), such as :group = To indicate an experimental design; full_mod = To fit a complete model or ref.batch = If we want to fit relative to a reference batch
    check_metadata(ribo,batch) # Verifies that the column specified by batch indeed exists in the ribo metadata and checks the consistency of the samples and their annotations
    matrix_ribo <- extract_data(ribo, "count", position_to_rownames = TRUE)
    # reorganize column according to metadata and convert DF to
    # matrix (otherwise, ComBat_seq won't work)
    matrix_ribo <- as.matrix(matrix_ribo[, c(ribo[["metadata"]][["samplename"]])])
    adjusted_matrix <- sva::ComBat_seq(matrix_ribo, batch = ribo[["metadata"]][[batch]], #applies CombatSeq, Indicates the batch variable from the metadata
                                       ...)
    ribo_updated <- .update_ribo_count_with_matrix(ribo, adjusted_matrix) # Replaces old counting data with adjusted data (adjusted_matrix)
    
    if (ribo_updated[["has_cscore"]]) { # if data already have computed c-score
      message("Recomputing c-score with the following parameters :",
              "\n- C-score method : ", ribo_updated[["cscore_method"]], 
              "\n- Flanking window : ", ribo_updated[["cscore_window"]],
              "\n")
      ribo_updated <- compute_cscore(ribo_updated, ribo_updated[["cscore_window"]], # then he recalculates it with the adjusted data
                                     ribo_updated[["cscore_method"]], ncores)

    }
    
    ribo_updated["combatSeq_count"] <- TRUE
    ribo_updated["col_used_combatSeq"] <- batch
    return(ribo_updated)

  
}

