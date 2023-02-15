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
#' show_rna_names(ribo = ribo_toy)
show_rna_names <- function(ribo = NULL) {
  if(!inherits(ribo,"RiboClass")) stop("A RiboClass object must be provided")
  RNA_names <- ribo[["rna_names"]][[2]]
  return(RNA_names)
}