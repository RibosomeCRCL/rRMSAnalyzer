#' Compute statistics (p-value for, Welch, Wilcoxon, Kruskal, ANOVA test) for C-score based on condition column
#'
#' @param test statistic test wanted student, anova or kruskal
#' @param condition_col column condition in metadata
#' @param ribo a RiboClass object
#'
#' @returns a dataframe
#' @export
#'
#' @examples 
#' data("ribo_toy")
#' res_pv(ribo = ribo_toy, test = "student", condition_col = "condition")
res_pv <- function(ribo = ribo, test = NULL, condition_col = NULL) {
  
  # Extract data
  data <- extract_data(ribo, only_annotated = TRUE, position_to_rownames = TRUE)
  
  # Transform data into long format excluding -annotated_sites
  data <- data %>%
    tibble::rownames_to_column("annotated_sites") %>%  # Add annotated_sites as column
    tidyr::pivot_longer(cols = -annotated_sites, names_to = "samplename", values_to = "c_score")  # keep annotated_sites outside pivoting
  
  # fusion with metadata
  data <- data %>%
    dplyr::left_join(ribo$metadata, by = "samplename") # left join between metadata and ribom_long
  
  #Erase lines with at least one NA
  datatest <- data %>% tidyr::drop_na({{condition_col}})
  
  # Determine test if not provided
  if (is.null(test)) {
    n_conditions <- length(unique(data[[condition_col]]))
    test <- if (n_conditions == 2) "student" 
    else if (n_conditions > 2) "anova"
    else message("Test can't be done with one group only")
  }
  
  # Test argument verification
  if (!test %in% c("student", "anova", "kruskal", "wilcoxon")) {
    stop("Test must be one of: 'student', 'anova', 'kruskal' or 'wilcoxon'")
  }
  
  #----------------------------------------------------------------------------
  #             Verify white space
  #----------------------------------------------------------------------------             
  if (any(grepl("\\s", ribo$metadata$samplename))) {
    stop(" Error samplename has whitespaces. Please rename your samplename.")
  }
  #----------------------------------------------------------------------------
  
  # Stats
  #create the table res_pv
  if (test %in% c("anova", "kruskal")) {
  res_pv <- data %>%
    dplyr::group_by(annotated_sites) %>%
    dplyr::summarise(
      p_value = {
        rep_table <- dplyr::cur_data()
        groups_table <- table(rep_table[[condition_col]])
        if (length(groups_table) >= 2 && all(groups_table >= 2)) {
          if (test == "anova") { # test == "anova"
            anova(lm(c_score ~ get(condition_col), data = rep_table))$"Pr(>F)"[1]
          } else { # test == "kruskal"
            kruskal.test(c_score ~ get(condition_col), data = rep_table)$p.value
          }
        } else {
          message(" Site ", unique(rep_table$annotated_sites),
                  " ignored : doesn't have at least 2 replicats.")
          NA_real_
        }
      }
    )
  #     p_value = anova(lm(c_score ~ get(condition_col)))$"Pr(>F)"[1])
  # } else {
  } else { #if not anova, nor kruskal  
    res_pv <- data %>%
      dplyr::group_by(annotated_sites) %>%
      #-------------------------------------------------------------------------
      dplyr::summarise(
        p_value = {
          # Check replicats
          rep_table <- dplyr::cur_data() # extract sub data
          groups_table <- table(rep_table[[condition_col]])
          # Verify number of group = 2 replicats number >= 2
          if (length (groups_table) == 2 && all(groups_table >= 2)) {
            if (test == "student") { # Welch here
              t.test(c_score ~ get(condition_col), var.equal = FALSE)$p.value #, var.equal = FALSE for Welch test, otherwise its student test
            } else if (test == "wilcoxon") {
                wilcox.test(c_score ~ get(condition_col))$p.value # Mann–Whitney U test because paired = FALSE (default)
              } else {
                  NA_real_
                }
        } else {
          message("Site ", unique(rep_table$annotated_sites), " ignored : doesn't have at least 2 replicats.")
          NA_real_
          }
        },
        delta_c_score = {
          rep_table <- dplyr::cur_data()
          groups_table <- table(rep_table[[condition_col]])
          if (length(groups_table) == 2 && all(groups_table >= 2)) {
            abs(mean(rep_table$c_score[rep_table[[condition_col]] == unique(rep_table[[condition_col]])[1]]) -
                  mean(rep_table$c_score[rep_table[[condition_col]] == unique(rep_table[[condition_col]])[2]]))
          } else {
            NA_real_
          }
        }
      )
}
      #-------------------------------------------------------------------------
            
  #           abs(mean(c_score[get(condition_col) == unique(get(condition_col))[1]]) -
  #   mean(c_score[get(condition_col) == unique(get(condition_col))[2]])) # delta C-score for student, wilcoxon and kruskal only
  # )
  # }
  
  # pval_adj for all tests
  res_pv$p_adj <- p.adjust(res_pv$p_value, method = "fdr")
  
  return(res_pv)
}
