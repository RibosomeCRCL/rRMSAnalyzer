#' Print Welch result in report
#'
#' @param ribo a RiboClass object
#' @param condition_col the condition column used for calculation in metadata
#' @param ctrl control in comparisons table
#' @param cases case in comparison table
#' @param pthr threshold of the p_value
#' @param cscore_cutoff threshold of c-score
#'
#' @returns a tag list
#' @export
#' @keywords internal
#' @examples 
#' data("ribo_toy")
#' data("human_methylated")
#' ribo_toy <- remove_ribo_samples(ribo_toy,c("RNA1", "RNA2"))
#' ribo_toy <- rename_rna(ribo_toy)  
#' ribo_toy <- annotate_site(ribo_toy,
#'                                 annot = human_methylated,
#'                                 anno_rna = "rRNA",
#'                                 anno_pos = "Position",
#'                                 anno_value = "Nomenclature")
#' comparisons <- tibble::tibble(comp = c("comp1"),
#'                               ctrl = c("cond1"),
#'                               cases = c("cond2") )  
#' ctrl <- comparisons$ctrl
#' cases <- comparisons$cases      
#' render_html_welch(ribo_toy, "condition", ctrl, cases, pthr = 0.05, cscore_cutoff = 0.05)
render_html_welch <- function(ribo, condition_col, ctrl, cases, pthr, cscore_cutoff) {
  ribo_filtered <- ribo
  kept_samples <- ribo_filtered$metadata %>%
    dplyr::filter(!!rlang::sym(condition_col) %in% c(cases, ctrl)) %>%
    dplyr::pull(samplename)
  
  ribo_filtered <- keep_ribo_samples(ribo_filtered, kept_samples)
  
  res_pv <- rRMSAnalyzer::res_pv(ribo = ribo_filtered, test = "student", condition_col = condition_col)
  
  res_pv_print <- res_pv
  if ("p_value" %in% names(res_pv_print)) {
    res_pv_print$p_value <- formatC(res_pv_print$p_value, format = "e", digits = 3)
  }
  if ("p_adj" %in% names(res_pv_print)) {
    res_pv_print$p_adj <- formatC(res_pv_print$p_adj, format = "e", digits = 3)
  }
  if ("delta_c_score" %in% names(res_pv_print)) {
    res_pv_print$delta_c_score <- round(res_pv_print$delta_c_score, 3)
  }
  
  # Creation of the table
  dt <- DT::datatable(res_pv_print,
                  rownames = FALSE,
                  extensions = c("Buttons"),
                  options = list(
                    pageLength = 20,
                    scrollX = TRUE,
                    dom = 'Bfrtip',
                    buttons = list(
                      list(extend = 'csv', filename = 'summary_table_of_Welch_test'),
                      list(extend = 'excel', filename = 'summary_table_of_Welch_test'),
                      list(extend = 'pdf', filename = 'summary_table_of_Welch_test')
                    )
                  ))

  w <- rRMSAnalyzer::plot_stat(ribo = ribo_filtered, site = NULL, res_pv = res_pv, pthr = pthr, condition_col = condition_col, cscore_cutoff = cscore_cutoff)
  empty_w <- list(plot_env = list(significant_sites = character(0)))
  # Conditionnal text
  text <- if (length(w[["plot_env"]][["significant_sites"]]) != 0) {
    glue::glue("This panel displays boxplots of C-scores for rRNA 2'Ome sites identified 
         as significantly different between the two biological conditions **{cases}** 
         and {ctrl}, based on Welch's t-test results. Significance was defined by 
         two criteria: adjusted p-value < 0.05 and |deltaC-score| > 0.05. Each boxplot 
         shows the distribution of C-scores per condition, based on the summary 
         table presented above. Colors are attributed to each biological group, 
         enabling visual comparison of 2'Ome levels between conditions.")
  } else {
    glue::glue("No significantly differentially rRNA 2'Ome sites were identified under 
         the current thresholds (Welch's t-test with adjusted p-value < 0.05 and 
         |deltaC-score| > 0.05).")
  }
  
  # Return all in tagList 
  output_list <- list(
    htmltools::tags$h4(glue::glue("{cases} vs {ctrl}")),
    dt,
    htmltools::tags$p(htmltools::HTML(text))
  )
  
  if (length(w[["plot_env"]][["significant_sites"]]) != 0) {
    # return plot in tagList and return NULL
    return(list(
      ui = htmltools::tagList(output_list),
      plot = w
    ))
  } else {
    return(list(
      ui = htmltools::tagList(output_list), plot = empty_w))
  }
}