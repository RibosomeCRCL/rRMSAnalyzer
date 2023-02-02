#' Keep only selected samples from a RiboClass Object.
#'
#' @param ribo A RiboClass object.
#' @param samples_to_keep A vector containing the names of samples to keep.
#'
#' @return A RiboClass with only selected samples in data and metadata.
#' @export
#'
#' @examples
#' data("ribo_toy")
#' ribo_test <- keep_ribo_samples(ribo = ribo_toy, samples_to_keep = c("S1", "S2"))

keep_ribo_samples <- function(ribo = NULL, samples_to_keep = NULL) {
  
  if(all(samples_to_keep %in% names(ribo[["data"]]))) {
    
    ribo[["data"]] <- ribo[["data"]][samples_to_keep]
    ribo[["metadata"]] <- ribo[["metadata"]][match(samples_to_keep, rownames(ribo[["metadata"]])),]
    
  } else {
    stop("Samples name should be from ", toString(names(ribo[["data"]])))
  }
  return(ribo)
}


#' Remove samples from a Riboclass object
#'
#' @param ribo A RiboClass object.
#' @param samples_to_delete A vector containing the names of samples to remove.
#'
#' @return A RiboClass without the removed samples in data and metadata.
#' @export
#'
#' @examples
#' data("ribo_toy")
#' ribo_test <- remove_ribo_samples(ribo = ribo_toy, samples_to_delete = c("S1", "S2"))

remove_ribo_samples <- function(ribo = NULL, samples_to_delete = NULL) {
  
  if(all(samples_to_delete %in% names(ribo[["data"]]))) {
    idx <- match(samples_to_delete, names(ribo[["data"]]))
    ribo[["data"]] <- ribo[["data"]][-idx]
    ribo[["metadata"]] <- ribo[["metadata"]][-match(samples_to_delete, rownames(ribo[["metadata"]])),]
    
  } else {
    stop("Samples name should be from ", toString(names(ribo[["data"]])))
  }
  return(ribo)
}
