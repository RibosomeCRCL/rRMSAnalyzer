#' Known 2'OMe positions in Humans' rRNA
#'
#' A dataset containing 112 2'OMe positions with associated data
#'
#' @format A data frame with 112 rows and 10 variables:
#' \describe{
#'   \item{Position}{nucleotide position on rRNA}
#'   \item{rRNA}{rRNA where this methylation is found}
#'   \item{Nomenclature}{name given to this modification}
#'   \item{NR_046235.Numbering}{Position on}
#'   \item{Sequence}{nucleotides sequence surround position}
#'   \item{SNORD}{ID of associated snoRNA}
#'   \item{Mode.of.coding}{how the SNORD is coded}
#'   \item{SNORD.host.gene}{SNORD's host gene}
#'   \item{Ensembl}{SNORD's Ensembl reference}
#'   \item{Nucleotide}{}
#'   ...
#' }
"human_methylated"

#' Suspected 2'OMe positions in Humans' rRNA
#'
#' A dataset containing 17 2'OMe positions with associated data
#'
#' @format A data frame with 17 rows and 10 variables:
#' \describe{
#'   \item{Position}{nucleotide position on rRNA}
#'   \item{rRNA}{rRNA where this methylation is found}
#'   \item{Nomenclature}{name given to this modification}
#'   \item{NR_046235.Numbering}{Position on}
#'   \item{Sequence}{nucleotides sequence surround position}
#'   \item{SNORD}{ID of associated snoRNA}
#'   \item{Mode.of.coding}{how the SNORD is coded}
#'   \item{SNORD.host.gene}{SNORD's host gene}
#'   \item{Ensembl}{SNORD's Ensembl reference}
#'   \item{Nucleotide}{}
#'   ...
#' }
"human_suspected"

#' RiboClass from a toy dataset
#'
#' A riboclass containing 10 samples + 2 reference RNA
"ribo_toy"