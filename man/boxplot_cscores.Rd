% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_boxplot.R
\name{boxplot_cscores}
\alias{boxplot_cscores}
\title{Plot boxplot representing the C-score values of all samples for each
individual annotated site.
Sites are sorted by their median.}
\usage{
boxplot_cscores(
  ribo,
  outlier = TRUE,
  sort_by = c("median", "iqr", "var")[1],
  horizontal = FALSE
)
}
\arguments{
\item{ribo}{A RiboClass object.}

\item{outlier}{Show boxplot outlier values.}

\item{sort_by}{Sort sites by median ("median", default) by variance ("var")
or IQR ("iqr").}

\item{horizontal}{Show boxplot horizontally.}
}
\value{
a ggplot geom_boxplot
}
\description{
Plot boxplot representing the C-score values of all samples for each
individual annotated site.
Sites are sorted by their median.
}
\examples{
data("ribo_toy")
data("human_methylated")
ribo_toy <- rename_rna(ribo_toy)
ribo_toy <- annotate_site(ribo_toy,human_methylated)
boxplot_cscores(ribo_toy)

}
