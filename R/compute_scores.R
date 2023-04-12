# (internal) Compute a c-score for all position of a single RNA
.compute_rna_cscore <- function(ds, flanking = 6, method) {
  # check that all parameters exist
  if (is.null(ds)) {
    stop("MISSING parameter. Please specify a data frame <ds>.")
  }
  
  # remove any cscore-related columns if they already exist
  ds[, c("flanking_median", "flanking_mad", "flanking_mean", "cscore")] <- list(NULL)
  # get the count column
  data_counts_col <- "count"
  
  ## c-score computation ##
  count <- ds[, data_counts_col]
  
  # First, we compute the window around each position (The latter
  # is excluded)
  if (method == "median") {
    flanking_values <- vapply(seq(1 + flanking, length(count)), function(x) {
      stats::median(count[c((x - flanking):(x - 1), (x + 1):(x + flanking))])
    }, numeric(1))
  } else if (method == "mean") {
    flanking_values <- vapply(seq(1 + flanking, length(count)), function(x) {
      mean(count[c((x - flanking):(x - 1), (x + 1):(x + flanking))])
    }, numeric(1))
  }
  # Because we started at the 1+flanking position, we append NA at
  # the beginning to match the original vector's size.
  flanking_values <- append(flanking_values, rep(NA, flanking), after = 0)
  
  if (method == "median") {
    ds[, "flanking_median"] <- flanking_values
    scorec_median_raw <- 1 - ds[, data_counts_col]/ds[, "flanking_median"]
    ds[, "cscore"] <- ifelse(scorec_median_raw < 0, 0, scorec_median_raw)
  } else if (method == "mean") {
    ds[, "flanking_mean"] <- flanking_values
    scorec_mean_raw <- 1 - ds[, data_counts_col]/ds[, "flanking_mean"]
    ds[, "cscore"] <- ifelse(scorec_mean_raw < 0, 0, scorec_mean_raw)
  }
  return(ds)
}

# (internal) Compute c-score for a single sample
.compute_sample_cscore <- function(sample_df = NULL, flanking = 6, method) {
  data_rna_col <- "rna"
  sample_df[, data_rna_col] <- as.factor(sample_df[, data_rna_col])
  RNA_counts_list <- split(sample_df, sample_df[, data_rna_col])
  sample_score <- lapply(RNA_counts_list, .compute_rna_cscore, flanking,
                         method)
  
  return(dplyr::bind_rows(sample_score))
  
}

#' Compute c-score for all samples in a riboclass
#' 
#' @description
#' 
# The C-score corresponds to the 2'Ome level at a RNA position. The
# C-score represents a drop in the end read coverage at a given
# position compared to the environmental coverage, as described by
# @birkedal2014. The C-score can be of 0 (i.e., no RNA molecule is
# 2'Ome at the position of interest), of 1 (i.e., all the RNA
# molecules are 2'Ome at the position of interest) and of ]0:1[
# (i.e., a mix of un-methylated and methylated RNA molecules).
#' 
#' In this package, the C-score is calculated for every position. As it can be
#' useful to find positions not yet identified as methylated.
#' 
#' This function is called automatically when loading data or after adjusting
#' biases in end read count data with ComBat-seq. It can be still called
#' manually to modify how the C-score is currently computed.
#' 
#' For each RNA, the first and last positions cannot be calculated if the local
#' coverage is shorter than the flanking argument. Their value will be NA instead.
#' 
#' @references Birkedal, U., Christensen-Dalsgaard, M., Krogh, N., Sabarinathan,
#' R., Gorodkin, J. and Nielsen, H. (2015),
#' Profiling of Ribose Methylations in RNA by High-Throughput Sequencing.
#' Angew. Chem. Int. Ed., 54: 451-455. https://doi.org/10.1002/anie.201408362
#' 
#'
#' @md
#' @param ribo A RiboClass object.
#' @param flanking Size of the local coverage.
#' @param method Computation method of the local coverage. Either 'median' or 'mean'.
#' @param ncores Number of cores to use in case of multithreading.
#' @return A RiboClass with c-score columns appended to each sample's data.
#' @export
#'
#' @examples
#' data('ribo_toy')
#' ribo_subsetted <- keep_ribo_samples(ribo_toy,'S1')
#' ribo_with_cscore_med <- compute_cscore(ribo_subsetted, ncores = 2)
#' ribo_with_cscore_mean <- compute_cscore(ribo_subsetted, method = 'mean', ncores = 2)
#' 
compute_cscore <- function(ribo = NULL, flanking = 6, method = "median",
                           ncores = 1) {
  if (!(method %in% c("median", "mean")))
    stop("method can be either \"mean\" or \"median\"")
  dt <- ribo["data"]  #we only need the counts to compute the score
  
  
  # Experimental : Multithreading is 3x faster than single-thread
  # TODO : implement
  if (ncores > 1) {
    samples_czscore <- parallel::mclapply(dt[["data"]], .compute_sample_cscore,
                                          flanking, method, mc.cores = ncores)
  } else {
    samples_czscore <- lapply(dt[["data"]], .compute_sample_cscore, flanking, method)
  }
  
  ribo[["data"]] <- samples_czscore
  ribo["cscore_window"] <- flanking
  ribo["cscore_method"] <- method
  ribo["has_cscore"] <- TRUE
  return(ribo)
}