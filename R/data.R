#' Known 2'Ome positions in Humans' rRNA
#'
#' A dataset containing 112 2'Ome positions with associated data
#' 
#' @usage data(human_methylated)
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
#'   \item{Nucleotide}{the nucleotide present at the position (A,T,G or C)}
#'   ...
#' }
"human_methylated"

#' Suspected 2'Ome positions in Humans' rRNA
#'
#' A dataset containing 17 2'Ome positions with associated data
#' @usage data(human_suspected)
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
#'   \item{Nucleotide}{the nucleotide present at the position (A,T,G or C)}
#'   ...
#' }
"human_suspected"

#' RiboClass from a toy dataset
#'
#' A riboclass containing 10 samples + 2 reference RNA.
#' 
#' Samples are from 4 different biological conditions ("condition" column in metadata).
#' The sequencing has been done in two different batches ("run" column in metadata).
#' Both batches have the same reference RNA, to detect technical bias.
#' @usage data(ribo_toy)
#' 
#' @format a RiboClass (S3 Class) with the following element
#' \describe{
#'   \item{data}{list of sample dataframe}
#'   \item{metadata}{metadata dataframe of all samples}
#'   \item{rna_names}{dataframe containing both original and modified RNAs name}
#'   }
"ribo_toy"