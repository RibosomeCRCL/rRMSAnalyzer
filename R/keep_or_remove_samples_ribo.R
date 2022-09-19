#' Keep only a selected number of samples
#'
#' @param ribo a RiboClass object, see constructor : 
#' \code{\link{create_riboclass}}
#' @param samples_to_keep A vector with the samples to keep
#'
#' @return a riboClass with only selected samples in data and metadata
#' @export
#'
#' @examples
#' ribo.test <- keep_ribo_samples(ribo = ribo_toy, samples_to_keep = c("4283", "3307"))

keep_ribo_samples <- function(ribo = NULL, samples_to_keep = NULL) {
  
  if(all(samples_to_keep %in% names(ribo[["data"]]))) {
    
    ribo[["data"]] <- ribo[["data"]][samples_to_keep]
    ribo[["metadata"]] <- ribo[["metadata"]][match(samples_to_keep, rownames(ribo[["metadata"]])),]
    
  } else {
    stop(paste("Samples name should be from", toString(names(ribo[["data"]]))))
  }
  return(ribo)
}


#' Remove samples from Riboclass object
#'
#' @param ribo a RiboClass
#' @param samples_to_delete A vector with the samples to remove
#'
#' @return a riboClass without the removed samples in data and metadata
#' @export
#'
#' @examples
#' ribo.test <- remove_ribo_samples(ribo = ribo_toy, samples_to_delete = c("4283", "3307"))

remove_ribo_samples <- function(ribo = NULL, samples_to_delete = NULL) {
  
  if(all(samples_to_delete %in% names(ribo[["data"]]))) {
    idx <- match(samples_to_delete, names(ribo[["data"]]))
    ribo[["data"]] <- ribo[["data"]][-idx]
    ribo[["metadata"]] <- ribo[["metadata"]][-match(samples_to_delete, rownames(ribo[["metadata"]])),]
    
  } else {
    stop(paste("Samples name should be from", toString(names(ribo[["data"]]))))
  }
  return(ribo)
}
