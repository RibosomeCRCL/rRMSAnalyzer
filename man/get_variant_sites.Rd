% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_variant_sites.R
\name{get_variant_sites}
\alias{get_variant_sites}
\title{Return the most/less variant sites of a dataframe or a RiboClass}
\usage{
get_variant_sites(
  ribo,
  n = NULL,
  type_of_variant = "most",
  only_annotated = TRUE
)
}
\arguments{
\item{ribo}{a RiboClass object.}

\item{n}{Number of top sites to return.}

\item{type_of_variant}{Either "less" or "most", to select the top n less
variant sites or the top n most variant sites respectively.}

\item{only_annotated}{(RiboClass only) Check variability only among
annotated sites. Ignored when df is a dataframe.}
}
\value{
A dataframe with the n most/less variant sites
}
\description{
Return the most/less variant sites of a dataframe or a RiboClass
}
\examples{
data("ribo_toy")
data("human_methylated")
ribo_toy <- rename_rna(ribo_toy)
ribo_toy <- annotate_site(ribo_toy,human_methylated)
get_variant_sites(ribo = ribo_toy)
}
