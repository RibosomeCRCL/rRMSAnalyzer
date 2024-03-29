% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_counts_fraction.R
\name{plot_counts_fraction}
\alias{plot_counts_fraction}
\title{Plot count distribution among RNAs for each sample}
\usage{
plot_counts_fraction(ribo)
}
\arguments{
\item{ribo}{A RiboClass object.}
}
\value{
a ggplot object
}
\description{
Plot count distribution among RNAs for each sample
}
\examples{
data("ribo_toy")
plot_counts_fraction(ribo_toy)
}
