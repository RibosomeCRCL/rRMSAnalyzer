#' Display RNA names
#'
#' @param ribo a RiboClass object
#'
#' @return
#' A vector with actual RNA names
#' @export
#'
#' @examples
#' data("ribo_toy")
#' show_RNA_names(ribo = ribo_toy)
show_RNA_names <- function(ribo = NULL) {
  if(is.null(ribo)) {stop("A ribo class object should be provided")}
  RNA_names <- ribo[["rna_names"]][[2]]
  return(RNA_names)
}