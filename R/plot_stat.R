#' Plot statistic of significant differential sites between experimental condition returned by the 
#' res.pv function 
#' @param ribo a RiboClass
#' @param site sites to plot
#' @param res_pv df of p values extracted from res_pv.R
#' @param pthr p value threshold
#' @param condition_col the condition column passed in params of the Rmd
#' @param p_cutoff cutoff of p-value under which sites are not taken in account
#' @param cscore_cutoff cutoff of c-score under which sites are not taken in account
#' @param adjust_pvalues_method Method to adjust p-value
#' @importFrom stats anova cor kruskal.test lm mad median p.adjust sd setNames t.test wilcox.test
#' @importFrom utils head
#' @export
#' @return a ggplot object
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
#' res_pv <- res_pv(ribo = ribo_toy, test = "kruskal", condition_col = "condition") 
#' plot_stat(ribo = ribo_toy, site = NULL, res_pv = res_pv, pthr = 0.01, 
#' condition_col = "condition", cscore_cutoff = 0.5)
plot_stat <- function(ribo = ribo,
                      p_cutoff = 0.05,
                      cscore_cutoff = 0.5,
                      adjust_pvalues_method = "fdr",
                      site = NULL, 
                      res_pv = res_pv, 
                      pthr = NULL, 
                      condition_col = NULL) {
  
  # Try to get params$condition_col if condition_col is NULL
  if (is.null(condition_col)) {
    if (exists("params") && !is.null(params$condition_col)) {
      condition_col <- params$condition_col
    } else {
      stop("Erreur : condition_col is not defined in parameters of the function or in the Rmd")
    }
  }
  
  #-------------------------------Case where plot can't be draw-----------------
  if (is.null(site) && (is.null(res_pv) || is.null(pthr))) {
    stop("Error : You need to specify at least one site")
  }
  if (!is.null(site) && is.null(res_pv) && !is.null(pthr)) {
    stop("Error : res_pv is required if pthr is given")
  }

  # ------------------------------Filter on sites-------------------------------
  # if site specified we use it without applying any threshold
  if (!is.null(site)) {
    significant_sites <- site
  } else if ("delta_c_score" %in% colnames(res_pv) ) { # if sites not specified, compute significant sites for welch and wilcoxon
    significant_sites <- res_pv %>% 
      dplyr::filter(p_adj < pthr, abs(delta_c_score - 1) >= cscore_cutoff) %>% 
      dplyr::pull(annotated_sites)
  } else {  # if sites not specified, compute significant sites for anova and kruskal
    warning("filter on p_adj")
    significant_sites <- res_pv %>% 
      dplyr::filter(p_adj < pthr) %>% 
      dplyr::pull(annotated_sites)
  }
  
  if (length(significant_sites) == 0) {
    return(ggplot() + 
             annotate("text", x = 4, y = 25, size = 8,
                      label = "No differential site found !") + 
             theme_void())
  }
  
  # ------------------------------Data formatting ------------------------------
  # Extract
  data <- extract_data(ribo, only_annotated = TRUE, position_to_rownames = TRUE)
  
  # Transform data into long format excluding -annotated_sites
  data <- data %>%
    tibble::rownames_to_column("annotated_sites") %>%  # Add annotated_sites as column
    tidyr::pivot_longer(cols = -annotated_sites, names_to = "samplename", values_to = "c_score")  # keep annotated_sites outside pivoting
  
  # Verify that specified sites exist
  sites_in_data <- unique(data$annotated_sites) #existing sites
  if (!all(significant_sites %in% sites_in_data)) { # if sites doesn't exist
    warning("given site in 'site = ...' doesn't exist in data. Please, check for typos in ", 
            paste(setdiff(significant_sites, sites_in_data), collapse = ", "))
  }
  
  # Left join between metadata and ribom_long
  data <- data %>%
    dplyr::left_join(ribo$metadata, by = "samplename") 
  

  # Add pval in data if res_pv is given 
  if (!is.null(res_pv)) { # keep only sites for which pval_adj is < pthr
    data <- data %>%
      dplyr::left_join(res_pv %>% dplyr::select(annotated_sites, p_value, p_adj), by = "annotated_sites")
  }
  
  # Keep only specified sites if necessary
  if(!is.null(site)) {
  data <- data %>% dplyr::filter(annotated_sites %in% significant_sites)
  }
  
  # -----------------------------Generate plot---------------------------------
  stat <- ggplot(data, aes(x = !!sym(condition_col), y = c_score, fill = !!sym(condition_col))) +
    geom_boxplot() +
    facet_wrap(~ annotated_sites) +
        theme_bw() +
    labs(title = "Statistical Analysis", x = "Condition", y = "C_score")
  
  # Add p_values on plot if res_pv given
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
