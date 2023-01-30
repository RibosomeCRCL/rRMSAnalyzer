#' print() to display basic informations about a RiboClass
#'
#' @param x a RiboClass
#' @param ... base print params
#'
#' @export
#' @keywords internal
print.RiboClass <- function(x,...){ 
  cat(as.character(paste("a RiboClass with", length(x[["data"]]),"samples and", nrow(x[["rna_names"]]), "RNA(s) :\n")))
   rna <- table(as.factor(x[["data"]][[1]][["rna"]]))
   rna_df <- as.data.frame(rna)
   rna_df <- rna_df[order(rna_df[,2]),]
   cat(paste0("Name : ",rna_df[,1] , ", length : ", rna_df[,2],"\n"))
}