#' Generate html report about batch effect in the working directory
#' 
#' @param ribo a riboclass
#' @param library_col library column in the metadata
#' @param project_name Name of the project
#' @export
report_qc <- function(ribo, library_col, project_name = "", output_dir = getwd()) {
    ribo_name <- deparse(substitute(ribo))
    ribo_to_check <- ribo #used in Rmarkdown template # nolint
    path <- system.file("rmd", package = "Riboscore")
    rmarkdown::render(paste0(path, "/check_bias.Rmd"),
                      params = list(
                                    library_col = library_col,
                                    ribo_name = ribo_name,
                                    project_name = project_name),
                      output_dir = output_dir)

}
