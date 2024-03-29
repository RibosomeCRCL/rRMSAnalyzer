% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/format_to_plot.R
\name{format_to_plot}
\alias{format_to_plot}
\title{Transform RiboClass data to a "ggplot-friendly" data frame.}
\usage{
format_to_plot(ribo, metadata_col = NULL, only_annotated = FALSE)
}
\arguments{
\item{ribo}{A RiboClass object.}

\item{metadata_col}{Metadata columns to add. Must be in RiboClass' metadata.}

\item{only_annotated}{Keep only sites that have been annotated.}
}
\value{
A ggplot-friendly dataframe with the following columns : 1) site; 2) sample; 3) cscore; ...) any metadata added
}
\description{
Turn a RiboClass into a dataframe that should be ready to be used in ggplot.
Metadata columns can be added if extra informations are needed.
}
