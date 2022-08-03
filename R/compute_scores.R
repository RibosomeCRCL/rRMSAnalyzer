#------------------------------------------

#' compute both C-score and Z-score for one RNA of one sample.
#'
#' @param ds dataframe of a given RNA for a given sample
#' @param flanking the number of flanking position to use for the window
#' @param data.counts.col column number where counts value are stored
#'
#' @return
#'
#' @examples
.compute_rna_cscore <- function(ds, flanking=6, method) {
  # check that all parameters exist
  if (is.null(ds)) {stop("MISSING parameter. Please specify a data frame <ds>.")}
  
  # remove any cscore-related columns if they already exist
  ds[,c("flanking_median","flanking_mad","flanking_mean","cscore")] <- list(NULL)
  # get the count column
   data_counts_col <- "count"
  
  # compute scores
  for (i in (flanking + 1) : (dim(ds)[1] - flanking)){
    #TODO : selection mean/median
    
  if(method == "median") {
    ds[i, "flanking_median"] <- stats::median(ds[c((i-flanking) : (i-1), (i+1) : (i + flanking)), data_counts_col]) # median
   # ds[i, "flanking_mad"] <- mad(ds[c((i-flanking) : (i-1), (i+1) : (i + flanking)), data_counts_col]) # mad
   # TODO : remove flanking mad ? 
    scorec_median_raw <- 1 - ds[, data_counts_col]/ds[, "flanking_median"]
    ds[, "cscore"] <- ifelse(scorec_median_raw < 0, 0, scorec_median_raw)
  }
    
 else if(method == "mean") {
    ds[i, "flanking_mean"] <- mean(ds[c((i-flanking) : (i-1), (i+1) : (i + flanking)), data_counts_col]) # mean
    scorec_mean_raw <- 1 - ds[, data_counts_col]/ds[, "flanking_mean"]
    ds[, "cscore"] <- ifelse(scorec_mean_raw < 0, 0, scorec_mean_raw)
 }
  }
  
  return(ds)
}


#' Title
#'
#' @param sample_df 
#' @param data_rna_col 
#' @param flanking 
#' @param data_counts_col 
#' @return
#'
#' @examples
.compute_sample_cscore <- function(sample_df=NULL,
                                      flanking=6,
                                    method) 
                                                   {
  data_rna_col = "rna"
  sample_df[,data_rna_col] <- as.factor(sample_df[,data_rna_col])
  RNA_counts_list <- split(sample_df, sample_df[,data_rna_col])
  sample_score <- lapply(RNA_counts_list, .compute_rna_cscore, flanking , method)
  
  return(dplyr::bind_rows(sample_score))
  
}

#' compute c-score for a given riboclass
#'
#' @param ribo a riboclass object 
#' @param flanking the window size around the position (the latter is excluded)
#' @param data_counts_col Name or position of the column containing count values
#' @param data_rna_col Name or position of the column containing the RNA
#' @param core number of core to use in case of multithreading
#' @return a riboclass with c-score columns added in data
#' @export
#'
#' @examples
#' ribo_with_cscore <- compute_cscore(ribo_toy)
#' 
#' @details 
#'\if{html}{\figure{cscore.png}{c-score visualised}}
#'\if{latex}{\figure{cscore.png}{options: width=0.5in}}' 
#' 
#' For each nucleotide position, a c-score is calculated by 
#' 
#' 
compute_cscore <- function(ribo=NULL, flanking=6,
                            method = "median",
                            use_multithreads = FALSE,
                            core = 8
                            ) {
  
   # Check if a cscore calculation has been already done on the ribo
   if( ribo["has_cscore"] == TRUE) {
     print("the riboClass has already a computed cscore")
   }

    dt = ribo["data"] #we only need the counts to compute the score

    
   # Experimental : Multithreading is 3x faster than single-thread
   # TODO : implement 
   if(use_multithreads) samples_czscore <- parallel::mclapply(dt[["data"]], .compute_sample_cscore, flanking, method,mc.cores = core)
   else samples_czscore <- lapply(dt[["data"]], .compute_sample_cscore, flanking, method)

     ribo[["data"]] <- samples_czscore
     ribo["cscore_window"] <- flanking
     ribo["cscore_method"] <- method
     ribo["has_cscore"] <- TRUE


  return(ribo)
  
  
}