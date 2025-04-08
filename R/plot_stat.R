#' Plot statistic
#'
#' @param ribo a RiboClass
#' @param site sites to plot
#' @param res_pv df of p values extracted from res_pv.R
#' @param pthr p value threshold
#' @param condition_col the condition column passed in params of the Rmd
#' @param p_cutoff cutoff of c-score under which sites are not taken in account
#' @param cscore_cutoff 
#' @param adjust_pvalues_method 
#'
#' @returns a boxplot
#' @export
#'
#' @examples plot_stat(data = ribo, site = NULL, res_pv = res_pv, pthr = 0.05)
plot_stat <- function(ribo = ribo_adj_annot,
                      p_cutoff = 1e-02,
                      cscore_cutoff = 0.5,
                      adjust_pvalues_method = "fdr",
                      site = NULL, 
                      res_pv = res_pv, 
                      pthr = NULL, 
                      condition_col = NULL) {
  library(ggplot2)
  
  # Try to get params$condition_col if condition_col is NULL
  if (is.null(condition_col)) {
    if (exists("params") && !is.null(params$condition_col)) {
      condition_col <- params$condition_col
    } else {
      stop("Erreur : condition_col is not defined in parameters of the function or in the Rmd")
    }
  }
  
    # Case where plot can't be draw
  if (is.null(site) && (is.null(res_pv) || is.null(pthr))) {
    stop("Error : You need to specify at least one site")
  }
  if (!is.null(site) && is.null(res_pv) && !is.null(pthr)) {
    stop("Error : res_pv is required if pthr is given")
  }

  # filter sites
  if (!is.null(res_pv) && !is.null(pthr)) {
    significant_sites <- res_pv %>% 
      dplyr::filter(p_adj < pthr, abs(fold_change - 1) >= cscore_cutoff) %>% 
      dplyr::pull(annotated_sites)
    if (length(significant_sites) == 0) {
      return(ggplot() + 
               annotate("text", x = 4, y = 25, size=8,
                        label = "No differential site found !") + 
               theme_void())
      }
  } else {
    significant_sites <- site
  }

  # Extract data
  data <- extract_data(ribo, only_annotated = TRUE, position_to_rownames = TRUE)
  
  # Transform data into long format excluding -annotated_sites
  data <- data %>%
    tibble::rownames_to_column("annotated_sites") %>%  # Add annotated_sites as column
    tidyr::pivot_longer(cols = -annotated_sites, names_to = "samplename", values_to = "c_score")  # keep annotated_sites outside pivoting
  
  data <- data %>%
    dplyr::left_join(ribo$metadata, by = "samplename") # left join between metadata and ribom_long
  
  data <- data %>%
    dplyr::filter(annotated_sites %in% significant_sites)
  
  # add pval if res_pv is given 
  if (!is.null(res_pv)) { # keep only sites for which pval_adj is < pthr
    data <- data %>%
      dplyr::left_join(res_pv %>% dplyr::select(annotated_sites, p_value, p_adj), by = "annotated_sites")
  }
  
  # Generate plot
  stat <- ggplot(data, aes(x = !!sym(condition_col), y = c_score, fill = !!sym(condition_col))) +
    geom_boxplot() +
    facet_wrap(~ annotated_sites) +
        theme_bw() +
    labs(title = "Statistical Analysis", x = "Condition", y = "C_score")
  
  # Add p_values if res_pv given
  if (!is.null(res_pv)) {
    stat <- stat + 
      geom_text(
        aes(x = 1.5, y = 0.3, 
            label = paste0("p=", signif(p_value, 3), "\n p_adj=", signif(p_adj, 3))),
        data = res_pv %>% dplyr::filter(annotated_sites %in% significant_sites),
        inherit.aes = FALSE,
        size = 3
      )
  }
  
return(stat)
} 
