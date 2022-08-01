#' Keep only a selected number of samples
#'
#' @param ribo 
#' @param samples_to_keep A vector with the samples to keep
#'
#' @return
#' @export
#'
#' @examples
#' ribo.test <- keep_ribo_samples(ribo = ribo, samples_to_keep = c("SAMPLE1", "SAMPLE2"))

keep_ribo_samples <- function(ribo = NULL, samples_to_keep = NULL) {
  
  if(all(samples_to_keep %in% names(ribo[["counts"]]))) {
    
    ribo[["counts"]] <- ribo[["counts"]][samples_to_keep]
    ribo[["metadata"]] <- ribo[["metadata"]][match(samples_to_keep, rownames(ribo[["metadata"]])),]
    
  } else {
    stop(paste("Samples name should be from", toString(names(ribo[["counts"]]))))
  }
  return(ribo)
}


#' Remove samples from Riboclass object
#'
#' @param ribo 
#' @param samples_to_delete A vector with the samples to keep
#'
#' @return
#' @export
#'
#' @examples
#' ribo.test <- remove_ribo_samples(ribo = ribo, samples_to_keep = c("SAMPLE1", "SAMPLE2"))

remove_ribo_samples <- function(ribo = NULL, samples_to_delete = NULL) {
  
  if(all(samples_to_delete %in% names(ribo[["counts"]]))) {
    idx <- match(samples_to_delete, names(ribo[["counts"]]))
    ribo[["counts"]] <- ribo[["counts"]][-idx]
    ribo[["metadata"]] <- ribo[["metadata"]][-match(samples_to_delete, rownames(ribo[["metadata"]])),]
    
  } else {
    stop(paste("Samples name should be from", toString(names(ribo[["counts"]]))))
  }
  return(ribo)
}
