#' Generate html report about possible batch effect for a given RiboClass.
#' 
#' @param ribo A RiboClass object.
#' @param library_col Library/run column in the metadata.
#' @param project_name Name of the project.
#' @param output_dir Path to output dir. Working directory is selected by default.
#' @export
report_qc <- function(ribo, library_col, project_name = "", output_dir = getwd()) {
    ribo_name <- deparse(substitute(ribo))
    ribo_to_check <- ribo #used in Rmarkdown template # nolint
    path <- system.file("rmd", package = "rRMSAnalyzer")
    rmarkdown::render(file.path(path, "check_bias.Rmd"),
                      params = list(
                                    library_col = library_col,
                                    ribo_name = ribo_name,
                                    project_name = project_name),
                      output_dir = output_dir)

}
