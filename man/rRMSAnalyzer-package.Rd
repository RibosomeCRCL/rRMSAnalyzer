% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rRMSAnalyzer.R
\docType{package}
\name{rRMSAnalyzer-package}
\alias{rRMSAnalyzer}
\alias{rRMSAnalyzer-package}
\title{rRMSAnalyzer package dedicated to analysis of RiboMethSeq data}
\description{
The rRMSAnalyzer package is designed to give easy-to-use functions
dedicated to RiboMethSeq (RMS) analyses. All functions revolve
around the RiboClass, a S3 class containing both RMS data and samples 
metadata.
}
\details{
The main functions are (in a standard workflow chronological order):

\itemize{
\item \code{\link{load_ribodata}} - import data and create a RiboClass. Objects of this class are used by all other functions
\item \code{\link{compute_cscore}} - compute C-score for all genomic positions in each sample
\item \code{\link{adjust_bias}} - adjust technical biases of your RiboMethSeq data with ComBat-seq
\item \code{\link{extract_data}} - get a dataframe of C-score for all samples and all positions in the RiboClass
\item Plots, e.g.: \code{\link{plot_pca}}; \code{\link{plot_coa}} or \code{\link{plot_heatmap}}.
}

To know how to start using this package, have a look at the vignette at 
\code{vignette("rRMSAnalyzer")}. It will give an example for most commands
and a workflow.


The code can be viewed at the GitHub repository :

\url{https://github.com/RibosomeCRCL/rRMSAnalyzer}
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/RibosomeCRCL/rRMSAnalyzer}
  \item \url{https://ribosomecrcl.github.io/rRMSAnalyzer/}
  \item Report bugs at \url{https://github.com/RibosomeCRCL/rRMSAnalyzer/issues}
}

}
\author{
Théo COMBE, Hermes PARAQINDES, Allyson MOUREAUX, Janice KIELBASSA
}
\keyword{package}
