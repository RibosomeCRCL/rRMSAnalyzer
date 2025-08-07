#' Function to generate HTML report using Rmarkdown
#' @description 
#' Note : if you want to create your own report template,
#' please make sure it has the following header : 
#' 
#' ---
#' title: "RMSAnalyzer : title of the report"
#' output: html_document
#' 
#' params:
#'  library_col: "library"
#'  condition_col: "condition"
#'  ribo_name: "ribo"
#'  project_name: "Unnamed project"
#' ---
#' 
#' @param ribo a RiboClass object.
#' @param rmdfile filename of the Rmarkdown template to use.
#' @param condition_col name or index of the column __in metadata__ containing the condition
#' @param library_col name or index of the column __in metadata__ containing the library
#' @param project_name Name of the project to display on the report
#' @param output_dir Path of output directory
#' @param ribo_name 
#' @param comments Path to a text file containing comments 
#' @return an html file
#' @keywords internal
report <- function(ribo, ribo_name,rmdfile, condition_col, library_col, project_name, output_dir, comments = NULL) {
  if (!is.null(comments) && file.exists(comments)) {
    comments_text <- readLines(comments, warn = FALSE)
  } else {
    comments_text <- "No comments given"
  }
  
    ribo_to_check <- ribo #used in Rmarkdown template # nolint
    path <- system.file("rmd", package = "rRMSAnalyzer")
    rmarkdown::render(file.path(path, rmdfile),
                      params = list(
                                    library_col = library_col,
                                    condition_col = condition_col,
                                    ribo_name = ribo_name,
                                    project_name = project_name,
                                    comments = comments_text),
                                    output_dir = output_dir)

}

#' Generate html report about possible batch effect for a given RiboClass.
#' 
#' @param ribo A RiboClass object.
#' @param library_col Library/run column in the metadata.
#' @param project_name Name of the project.
#' @param comments Path to a text file containing comments
#' @param output_dir Path to output dir. Working directory is selected by default.
#' @return an html file
#' @export
#' @examples 
#' data('ribo_toy')
#' report_qc(ribo_toy, "run")
report_qc <- function(ribo, library_col, project_name = "Unnamed project", output_dir = getwd(), comments = NULL) {
  ribo_name <- deparse(substitute(ribo))
  
  report(ribo, ribo_name,"quality_control.Rmd",condition_col = NULL,
         library_col = library_col,
         project_name = project_name,
         output_dir = output_dir,
         comments = comments)
}
