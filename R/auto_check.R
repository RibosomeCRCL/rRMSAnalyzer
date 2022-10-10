#' Generate html report about batch effect in the working directory
#' 
#' @param ribo a riboclass
#' @param library_col library column in the metadata
#' @export
batch_report <- function(ribo, library_col) {
    ribo_to_check <- ribo #used in Rmarkdown template
    path <- system.file("rmd", package = "Riboscore")
    rmarkdown::render(paste0(path,"/check_bias.Rmd"),
                      params = list(library_col = library_col),
                      output_dir = getwd())
}