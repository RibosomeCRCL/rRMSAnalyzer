% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/keep_or_remove_samples_ribo.R
\name{keep_ribo_samples}
\alias{keep_ribo_samples}
\title{Keep only selected samples from a RiboClass Object.}
\usage{
keep_ribo_samples(ribo = NULL, samples_to_keep = NULL)
}
\arguments{
\item{ribo}{A RiboClass object.}

\item{samples_to_keep}{A vector containing the names of samples to keep.}
}
\value{
A RiboClass with only selected samples in data and metadata.
}
\description{
Keep only selected samples from a RiboClass Object.
}
\examples{
data("ribo_toy")
ribo_test <- keep_ribo_samples(ribo = ribo_toy, samples_to_keep = c("S1", "S2"))
}
