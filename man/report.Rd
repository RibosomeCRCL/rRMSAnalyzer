% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/report.R
\name{report}
\alias{report}
\title{Function to generate HTML report using Rmarkdown}
\usage{
report(
  ribo,
  ribo_name,
  rmdfile,
  condition_col,
  library_col,
  project_name,
  output_dir,
  comments = NULL
)
}
\arguments{
\item{ribo}{a RiboClass object.}

\item{ribo_name}{}

\item{rmdfile}{filename of the Rmarkdown template to use.}

\item{condition_col}{name or index of the column __in metadata__ containing the condition}

\item{library_col}{name or index of the column __in metadata__ containing the library}

\item{project_name}{Name of the project to display on the report}

\item{output_dir}{Path of output directory}

\item{comments}{Path to a text file containing comments}
}
\description{
Note : if you want to create your own report template,
please make sure it has the following header : 

---
title: "RMSAnalyzer : title of the report"
output: html_document

params:
 library_col: "library"
 condition_col: "condition"
 ribo_name: "ribo"
 project_name: "Unnamed project"
---
}
\keyword{internal}
