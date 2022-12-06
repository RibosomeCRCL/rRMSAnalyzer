#' rRMSAnalyzer package dedicated to analysis pf RiboMethSeq data
#' 
#' The rRMSAnalyzer package is designed to give easy-to-use functions
#' dedicated to RiboMethSeq (RMS) analyses. All functions revolve
#' around the riboClass, a S3 class containing both RMS data and 
#' metadata.
#'
#' The main functions are (in a standard workflow chronological order):
#'
#' \itemize{
#' \item \code{\link{create_riboclass}} - import data and create a riboClass. It is used by all other functions
#' \item \code{\link{compute_cscore}} - compute C-score for all genomic positions in your samples
#' \item \code{\link{adjust_bias}} - adjust technical biases of your RiboMethSeq data with combatSeq
#' \item \code{\link{extract_data}} - get a dataframe of c-score for all samples and all positions in your riboClass
#' \item Plots, e.g.: \code{\link{plot_PCA}}.
#' }
#' 
#' To know how to start using this package, you can have a look at the vignette at 
#' \code{vignette("analyses-with-rRMSAnalyzer")}. It will give an example for most commands
#' and a workflow.
#' 
#'
#' The code can be viewed at the GitHub repository :
#' 
#' \url{https://github.com/...}
#' 
#'
#'
#' @author Théo COMBE, Hermes PARAQINDES, Janice KIELBASSA
#' 
#' @docType package
#' @name rRMSAnalyzer-package
#' @aliases rRMSAnalyzer-package
#' @keywords package
NULL