#' Compute statistics for c-score based on condition column
#'
#' @param ribo a RiboClass object
#' @param test_name name of the test to perform
#' @param col_condition condition column on which perform the tests
#'
#' @returns a dataframe
#' @export
#'
#' @examples res_pv(ribo = ribo_adj_annot, test_name = t.test, condition_col = col)
res_pv <- function(ribo = ribo_adj_annot, test = "student", condition_col = PR3) {
  # Test arg verification
  if(!test %in% c("student", "anova", "kruskal")) {
    stop("Le test doit être 't.test', 'anova' ou 'kruskal.test'")
  }
  
  # Extract data
  data <- extract_data(ribo, only_annotated = TRUE, position_to_rownames = TRUE)
  
  # Transform data into long format excluding -annotated_sites
  data <- data %>%
    tibble::rownames_to_column("annotated_sites") %>%  # Add annotated_sites as column
    tidyr::pivot_longer(cols = -annotated_sites, names_to = "samplename", values_to = "c_score")  # keep annotated_sites outside pivoting
  
  # fusion with metadata
  data <- data %>%
    dplyr::left_join(ribo$metadata, by = "samplename") # left join between metadata and ribom_long
  
  # Stats
  res_pv <- data %>%
    dplyr::group_by(annotated_sites) %>%
    dplyr::summarise(
      p_value = ifelse(test == "student", t.test(c_score ~ get(condition_col))$p.value,
                       ifelse(test == "anova", anova(lm(c_score ~ get(condition_col)))$"Pr(>F)"[1],
                              kruskal.test(c_score ~ get(condition_col))$p.value)),
      fold_change = mean(c_score[get(condition_col) == unique(get(condition_col))[1]]) /
        mean(c_score[get(condition_col) == unique(get(condition_col))[2]])
    )
  
  # pval_adj
  res_pv$p_adj <- p.adjust(res_pv$p_value, method = "fdr")
  
  return(res_pv)
}
