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
#' @param specie studied specie of the project
#' @param ribo_name name of RiboClass
#' @param rmdfile filename of the Rmarkdown template to use.
#' @param condition_col name or index of the column __in metadata__ containing the condition
#' @param library_col name or index of the column __in metadata__ containing the library
#' @param project_name Name of the project to display on the report
#' @param output_dir Path of output directory
#' @param comments Path to a text file containing comments 
#'
#' @return an html file
#' @keywords internal
report <- function(ribo, specie = "not specified", ribo_name, rmdfile, condition_col, library_col, project_name, output_dir, comments = NULL) {
  if (!is.null(comments) && file.exists(comments)) {
    comments_text <- readLines(comments, warn = FALSE)
  } else {
    comments_text <- "No comments given"
  }
  if (!specie %in% c("human","murine")) {
    stop("please specify the studied specie as 'human' or 'murine'")
  }
    ribo_to_check <- ribo #used in Rmarkdown template # nolint
    path <- system.file("rmd", package = "rRMSAnalyzer")
   
    rmarkdown::render(file.path(path, rmdfile),
                      params = list(specie = specie,
                                    library_col = library_col,
                                    condition_col = condition_col,
                                    ribo_name = ribo_name,
                                    project_name = project_name,
                                    comments = comments_text),
                                    output_dir = output_dir)

}

#' Generate html report about possible batch effect for a given RiboClass.
#' @inheritParams report
#' @return an html file
#' @export
#' @examples 
#' data('ribo_toy')
#' report_qc(ribo_toy, "human", "run")
report_qc <- function(ribo, specie = "not specified", library_col, project_name = "Unnamed project", output_dir = getwd(), comments = NULL) {
  ribo_name <- deparse(substitute(ribo))
  
  report(ribo, specie = specie,
         ribo_name,"quality_control.Rmd",condition_col = NULL,
         library_col = library_col,
         project_name = project_name,
         output_dir = output_dir,
         comments = comments)
}


#' Generate html report about the analysis of 2’Ome profile
#' @inheritParams report
#' @export
#' @examples 
#' data('ribo_toy')
#' data("human_methylated")
#' ribo_adj_small <- remove_ribo_samples(ribo_toy, c("RNA1","RNA2"))
#' ribo_adj_small_name <- rename_rna(ribo_adj_small, 
#' new_names = c("5S", "5.8S", "18S", "28S"))
#' ribo_adj_annot <- annotate_site(ribo_adj_small_name, 
#'                                 annot = human_methylated,
#'                                 anno_rna = "rRNA",
#'                                 anno_pos = "Position",
#'                                 anno_value = "Nomenclature") 
#'report_2ome_sites(ribo = ribo_adj_annot, specie = "human", 
#'condition_col = "condition", project_name = "name", comments = "./comment_2ome_comp1.Rmd") 
report_2ome_sites <- function(ribo, specie = "not specified", condition_col, project_name = "Unnamed project", output_dir = getwd(), comments = NULL) {
  ribo_name <- deparse(substitute(ribo))
  
  # Before making a report, we must check if the ribo has been annotated.
  # plot_diff_sites only compares between annotated sites.
  
  if(all(is.na(ribo[["data"]][[1]][["site"]]))) {
    stop("Please annotate your RiboClass before using report_2ome_sites !")
  }
  
  report(ribo, 
         specie = specie, 
         ribo_name,
         "2ome_analysis.Rmd",
         condition_col = condition_col,
         library_col = NULL,
         project_name = project_name,
         output_dir = output_dir, 
         comments = comments)
}

#' Generate html report about the statistical analysis of difference in 2’Ome 
#' level at each site
#' @inheritParams report
#' @param comparisons data.frame or tibble with columns 'ctrl' and 'cases', 
#' defining the pairwise comparisons to run.
#' @importFrom magrittr %>%
#' @export
#' @examples
#' library(dplyr)
#' data('ribo_toy')
#' data("human_methylated")
#' ribo_adj_small <- remove_ribo_samples(ribo_toy, c("RNA1","RNA2"))
#' ribo_adj_small_name <- rename_rna(ribo_adj_small, 
#' new_names = c("5S", "5.8S", "18S", "28S"))
#' ribo_adj_annot <- annotate_site(ribo_adj_small_name, 
#'                                 annot = human_methylated,
#'                                 anno_rna = "rRNA",
#'                                anno_pos = "Position",
#'                                anno_value = "Nomenclature")
#' comparisons <- tibble::tibble(
#' comp = c("comp1"),
#' ctrl = c("cond1"),
#' cases = c("cond2"))
#' kept_samples <- ribo_adj_annot$metadata %>% 
#'   dplyr::filter(!is.na(comp1)) %>% 
#'   dplyr::pull(samplename)
#' ribo_adj_annot_comp1 <- keep_ribo_samples(ribo_adj_annot, kept_samples)
#' report_diff_sites(ribo = ribo_adj_annot_comp1, specie = "human", 
#' condition_col = "comp1", project_name = "name", comparisons = comparisons, 
#' comments = "./comment_diff_site.Rmd") 

report_diff_sites <- function(ribo, specie = "not specified", ribo_name, rmdfile, condition_col, library_col, project_name = "Unnamed project", output_dir = getwd(), comparisons = NULL, comments = NULL) {
  
  ribo_name <- deparse(substitute(ribo))
  
  if (!is.null(comments) && file.exists(comments)) {
    comments_text <- readLines(comments, warn = FALSE)
  } else {
    comments_text <- "No comments given"
  }
  
  ribo_to_check <- ribo #used in Rmarkdown template # nolint
  path <- system.file("rmd", package = "rRMSAnalyzer")
  
  
  # Before making a report, we must check if the ribo has been annotated.
  # plot_diff_sites only compares between annotated sites.
  if (all(is.na(ribo[["data"]][[1]][["site"]]))) {
    stop("Please annotate your RiboClass before using report_diff_sites !")
  }
  
  if (is.null(comparisons)) {
    stop("You must provide a 'comparisons' table with columns 'ctrl' and 'cases'")
  }
  
  params_list <- list(specie = specie,
                      library_col   = NULL,
                      condition_col = condition_col,
                      ribo_name     = ribo_name,
                      project_name  = project_name,
                      comments      = comments_text,
                      comparisons   = comparisons
  )
  
  rmarkdown::render(
    file.path(path, "diff_sites.Rmd"),
    params = params_list,
    output_dir = output_dir
  )
}
