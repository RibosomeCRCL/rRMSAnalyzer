#' compute both C-score and Z-score for one RNA of one sample.
#'
#' @param ds dataframe of a given RNA for a given sample
#' @param flanking the number of flanking position to use for the window
#' @param method either "median" or "mean", depending on how the cscore must be computed.
#' 
#' @return a single RNA dataframe with cscore and flanking columns
#' @keywords internal

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


#' compute cscore for a single sample
#'
#' @param sample_df the whole a dataframe for a single sample (all rRNA combined)
#' @param flanking the number of flanking position to use for the window
#' @param method either "median" or "mean", depending how the cscore must be computed.
#' @return a sample dataframe with cscore and flanking columns
#'
#' @keywords internal
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

#' compute c-score for all samples in a riboclass
#' 
#' @description
#' 
#' The C-score corresponds to the 2'Ome level at a rRNA position known to be methylated. The C-score represents a drop in end read coverage at a given position compared to the environmental coverage as described by Birkedal et al, 2015. This score can have a value between **0** (never methylated) and **1** (always methylated).
#' 
#' In this package, the C-score is calculated for every position. As it can be useful to find positions not yet identified as methylated.
#' 
#' 
#' For each RNA, The first and last positions cannot be calculated due to window's size. Their value will be NA instead.
#' 
#' @references Birkedal, U., Christensen-Dalsgaard, M., Krogh, N., Sabarinathan, R., Gorodkin, J. and Nielsen, H. (2015), Profiling of Ribose Methylations in RNA by High-Throughput Sequencing. Angew. Chem. Int. Ed., 54: 451-455. https://doi.org/10.1002/anie.201408362
#' 
#'
#' @md
#' @param ribo a riboclass object, see constructor : 
#' \code{\link{create_riboclass}}
#' @param flanking the window size around the position (the latter is excluded)
#' @param method either "median" or "mean", depending how the cscore must be computed.
#' @param ncores number of ncores to use in case of multithreading
#' @return a riboclass with c-score related columns appended to each sample's data.
#' @export
#'
#' @examples
#' data("ribo_toy")
#' ribo_subsetted <- keep_ribo_samples(ribo_toy,"4283")
#' ribo_with_cscore_med <- compute_cscore(ribo_subsetted, ncores = 2)
#' ribo_with_cscore_mean <- compute_cscore(ribo_subsetted, method = "mean", ncores = 2)
#' 
compute_cscore <- function(ribo=NULL, flanking=6,
                            method = "median",
                            ncores = 1
                            ) {
    if(!(method %in% c("median","mean"))) stop("method can be either \"mean\" or \"median\"")
    dt = ribo["data"] #we only need the counts to compute the score

    
   # Experimental : Multithreading is 3x faster than single-thread
   # TODO : implement 
   if(ncores > 1) samples_czscore <- parallel::mclapply(dt[["data"]], .compute_sample_cscore, flanking, method,mc.cores = ncores)
   else samples_czscore <- lapply(dt[["data"]], .compute_sample_cscore, flanking, method)

     ribo[["data"]] <- samples_czscore
     ribo["cscore_window"] <- flanking
     ribo["cscore_method"] <- method
     ribo["has_cscore"] <- TRUE


  return(ribo)
  
  
}