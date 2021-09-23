#------------------------------------------

#' Calculate both C-score and Z-score for one RNA of one sample.
#'
#' @param ds 
#' @param flanking 
#' @param data.counts.col 
#' @param data.position.col 
#'
#' @return
#' @export
#'
#' @examples
calculate_score_by_RNA <- function(ds, flanking=6,
                                   data_counts_col=4,
                                   data_position_col=1) {
  # check that all parameters exist
  if (is.null(ds)) {stop("MISSING parameter. Please specify a data frame <ds>.")}
  ds[, "mean"] <- NA
  ds[, "median"] <- NA
  ds[, "mad"] <- NA
  
  # compute scores
  for (i in (flanking + 1) : (dim(ds)[1] - flanking)){
    #TODO : selection mean/median
    ds[i, "mean"] <- mean(ds[c((i-flanking) : (i-1), (i+1) : (i + flanking)), data_counts_col]) # mean
    ds[i, "median"] <- median(ds[c((i-flanking) : (i-1), (i+1) : (i + flanking)), data_counts_col]) # median
    ds[i, "mad"] <- mad(ds[c((i-flanking) : (i-1), (i+1) : (i + flanking)), data_counts_col]) # mad
    
  }
  ds[, "dist2medInMad"] <- (ds[,"median"] - ds[, data_counts_col])/ds[,"mad"]
  scorec_median_raw <- 1 - ds[, data_counts_col]/ds[, "median"]
  ds[, "ScoreC.Median.net"] <- ifelse(scorec_median_raw < 0, 0, scorec_median_raw)
  
  return(ds)
}


#' Title
#'
#' @param sample_df 
#' @param data_rna_col 
#' @param flanking 
#' @param data_counts_col 
#' @param data_position_col 
#'
#' @return
#' @export
#'
#' @examples
calculate_score_by_sample <- function(sample_df=NULL,
                                      data_rna_col = 1,
                                      flanking=6, 
                                      data_counts_col=3,
                                      data_position_col=1) {
  
  sample_df[,data_rna_col] <- as.factor(sample_df[,data_rna_col])
  RNA_counts_list <- split(sample_df, sample_df[,data_rna_col])
  sample_score <- lapply(RNA_counts_list, calculate_score_by_RNA, flanking = flanking, data_counts_col = data_counts_col, data_position_col = data_position_col)
  
  return(dplyr::bind_rows(sample_score))
  
}

#' Calculate score for a whole sample cohort
#'
#' @param dt 
#' @param flanking 
#' @param data_counts_col 
#' @param data_position_col 
#'
#' @return
#' @export
#'
#' @examples
calculate_score <- function(dt=NULL, flanking=6,
                            data_rna_col = 1,
                            data_counts_col=4,
                            use_multithreads = F) {
  
  
  if(class(dt) == "RiboClass") {
    ribo = dt
    dt = dt["raw_counts"] #we only need the counts to calculate the score
    
  }
   # Experimental : Multithreading is 3x faster than single-thread
   # TODO : implement 
   if(use_multithreads) samples_czscore <- BiocParallel::bplapply(dt[["raw_counts"]], calculate_score_by_sample)
   else samples_czscore <- lapply(dt[["raw_counts"]], calculate_score_by_sample, data_counts_col = data_counts_col, data_rna_col = data_rna_col)
  
   if(exists("ribo")) {
     ribo[["raw_counts"]] <- samples_czscore
     return(ribo)
   }
  return(samples_czscore)
  
  
}