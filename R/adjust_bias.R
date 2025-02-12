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
#' 
adjust_bias <- function(ribo, batch,ncores = 1, ...) { # ... permet de passer des arguments supplémentaires à sva::ComBat_seq(), comme :group = Pour indiquer un design expérimental ; full_mod = Pour ajuster un modèle complet ou ref.batch = Si on veut ajuster par rapport à un batch de référence
    check_metadata(ribo,batch) # Vérifie que la colonne spécifiée par batch existe bien dans les métadonnées de ribo et vérifie la cohérence des échantillons et de leurs annotations
    matrix_ribo <- extract_data(ribo, "count", position_to_rownames = TRUE)
    # reorganize column according to metadata and convert DF to
    # matrix (otherwise, ComBat_seq won't work)
    matrix_ribo <- as.matrix(matrix_ribo[, c(ribo[["metadata"]][["samplename"]])])
    adjusted_matrix <- sva::ComBat_seq(matrix_ribo, batch = ribo[["metadata"]][[batch]], #applique CombatSeq, Indique la variable de batch depuis les métadonnées
                                       ...)
    ribo_updated <- .update_ribo_count_with_matrix(ribo, adjusted_matrix) # Remplace les anciennes données de comptage par les données ajustées (adjusted_matrix)
    
    if (ribo_updated[["has_cscore"]]) { # si les données ont déjà un c-score de calculé
      message("Recomputing c-score with the following parameters :",
              "\n- C-score method : ", ribo_updated[["cscore_method"]], 
              "\n- Flanking window : ", ribo_updated[["cscore_window"]],
              "\n")
      ribo_updated <- compute_cscore(ribo_updated, ribo_updated[["cscore_window"]], # alors il le recalcule avec les données ajustées
                                     ribo_updated[["cscore_method"]], ncores)

    }
    
    ribo_updated["combatSeq_count"] <- TRUE
    ribo_updated["col_used_combatSeq"] <- batch
    return(ribo_updated)

  
}

