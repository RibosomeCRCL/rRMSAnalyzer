#' Compute PCA from a c-score matrix
#'
#' @keywords internal
#'
#' @param cscore.matrix matrix of c-score extracted from a RiboClass with ' \code{\link{extract_data}}
#' @return dudi.pca object
#'
.calculate_pca <- function(cscore.matrix = NULL) {
  pca.res <- ade4::dudi.pca(t(cscore.matrix[stats::complete.cases(cscore.matrix),]), 
                            scannf = FALSE, 
                            nf = 5)
  
  return(pca.res)
}