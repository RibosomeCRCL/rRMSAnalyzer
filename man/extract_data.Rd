% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_data.R
\name{extract_data}
\alias{extract_data}
\title{Aggregate results into a single matrix}
\usage{
extract_data(
  ribo,
  col = "cscore",
  position_to_rownames = FALSE,
  only_annotated = FALSE
)
}
\arguments{
\item{ribo}{A RiboClass object.}

\item{col}{Column in data you want extract data from (cscore or count).}

\item{position_to_rownames}{If true, position will be included as a rowname. They will in a new column otherwise.}

\item{only_annotated}{If true, return a dataframe with only annotated sites. Return all sites otherwise.}
}
\value{
A dataframe.
}
\description{
For a given column in data, this function will generate a dataframe with all samples.
Exports all positions (if only_annotated is false) or only positions with a siteID (if only_annotated is true).
}
\examples{
data("ribo_toy")
count_df <- extract_data(ribo_toy,"count")
}
