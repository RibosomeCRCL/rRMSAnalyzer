#' Get a summary of statistics for differential 2'ome-sites analysis
#'
#' @param ribo a RiboClass
#' @param pthr p-value threshold
#' @param condition_col the condition column used for calculation in metadata
#' @param cscore_cutoff cutoff of c-score under which sites are not taken in account
#' @param comparisons table given by the user
#'
#' @returns a data frame
#' @export
#'
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
#' comparisons <- tibble::tibble(
#'                               comp = c("comp1"),
#'                               ctrl = c("cond1"),
#'                               cases = c("cond2") )                                
#' get_diff_sites_summary(ribo = ribo_toy, pthr = 0.05, condition_col = "condition",
#'  cscore_cutoff = 0.05, comparisons = comparisons)
                      
get_diff_sites_summary <- function(ribo = ribo, pthr = 0.05, condition_col = condition_col, cscore_cutoff = 0.05, comparisons = comparisons){
  
  #extract metadata
  qcdata <- ribo[[2]] 
  
  #----------------------------------------ANOVA test---------------------------
  if(length(unique(qcdata[[condition_col]])) > 2){
  compute_pval_anova <- rRMSAnalyzer::compute_pval(ribo = ribo, test = "anova", condition_col = condition_col) 
  a <-rRMSAnalyzer::plot_stat(ribo = ribo, compute_pval = compute_pval_anova, pthr = pthr, condition_col = condition_col, cscore_cutoff = cscore_cutoff)
  sites_anova <- a[["plot_env"]][["significant_sites"]]
  }
  
  #----------------------------------------Welch test---------------------------
  w_list <- list() 
  for(i in seq_len(nrow(comparisons))) {
    ctrl <- comparisons$ctrl[i]
    cases <- comparisons$cases[i]
    
    ribo_filtered <- ribo
    kept_samples <- ribo$metadata %>%
      dplyr::filter(!!sym(condition_col) %in% c(cases, ctrl))%>%
      dplyr::pull(samplename)
    
    ribo_filtered <- keep_ribo_samples(ribo_filtered,kept_samples)
    
    compute_pval_welch <- rRMSAnalyzer::compute_pval(ribo = ribo_filtered, test = "student", condition_col = condition_col) 
    
    w <- rRMSAnalyzer::plot_stat(ribo = ribo_filtered, compute_pval = compute_pval_welch, pthr = pthr, condition_col = condition_col, cscore_cutoff = cscore_cutoff)
    w_list[[i]] <- w 
  }
  
  #----------------------------------------Kruskal test-------------------------
  if(length(unique(qcdata[[condition_col]])) > 2){
    compute_pval_kruskal <- rRMSAnalyzer::compute_pval(ribo = ribo, test = "kruskal", condition_col = condition_col)
    k <-rRMSAnalyzer::plot_stat(ribo = ribo, compute_pval = compute_pval_kruskal, pthr = pthr, condition_col = condition_col, cscore_cutoff = cscore_cutoff)
    sites_kw <- k[["plot_env"]][["significant_sites"]]
  }
  
  #----------------------------------------Wilcoxon test------------------------
  x_list <- list()
  for(i in seq_len(nrow(comparisons))) {
    ctrl <- comparisons$ctrl[i]
    cases <- comparisons$cases[i]
    
    ribo_filtered_x <- ribo
    kept_samples_x <- ribo$metadata %>%
      dplyr::filter(!!sym(condition_col) %in% c(cases, ctrl))%>%
      dplyr::pull(samplename)
    
    ribo_filtered_x <- keep_ribo_samples(ribo_filtered_x,kept_samples_x)
    
    compute_pval_wilcoxon <- rRMSAnalyzer::compute_pval(ribo = ribo_filtered_x, test = "wilcoxon", condition_col = condition_col)
    
    x <- rRMSAnalyzer::plot_stat(ribo = ribo_filtered_x, compute_pval = compute_pval_wilcoxon, pthr = pthr, condition_col = condition_col, cscore_cutoff = cscore_cutoff)
    x_list[[i]] <- x 
  }
  
  #-------------------------------------Summary table creation------------------
  # Init list
  comparison_results <- list()
  
  # Fill each comparison
  for (i in seq_len(nrow(comparisons))) {
    ctrl <- comparisons$ctrl[i]
    case <- comparisons$cases[i]
    comp_id <- paste0(case, "_vs_", ctrl)
    
    # Extract significant sites from tests
    sites_welch <- unique(unlist(lapply(w_list, function(w) w[["plot_env"]][["significant_sites"]]))) 
    sites_wilcox <- unique(unlist(lapply(x_list, function(x) x[["plot_env"]][["significant_sites"]])))
    
    # Add to the list
    comparison_results[[comp_id]] <- list(
      Welch = sites_welch,
      Wilcoxon = sites_wilcox
    )
  }
  
  # Init list
  results <- list()
  
  # Identify all significant sites in at least on test
  if(length(unique(qcdata[[condition_col]])) > 2){
    sites_anova <- a[["plot_env"]][["significant_sites"]]
    sites_kw <- k[["plot_env"]][["significant_sites"]]
  }
  
  # all sites
  if(length(unique(qcdata[[condition_col]])) > 2){
  all_sites <- unique(c(sites_anova, sites_kw, sites_welch, sites_wilcox))
  } else {
    all_sites <- unique(c(sites_welch, sites_wilcox))
  }
  
  # first column
  results$Site <- all_sites
  
  # Add Anova et Kruskal column if groups > 2 
  if (length(unique(qcdata[[condition_col]])) > 2) {
    results$Multiple_comparison_Anova <- ifelse(all_sites %in% sites_anova, "Significant", "NS")
    results$Multiple_comparison_Kruskal_Wallis <- ifelse(all_sites %in% sites_kw, "Significant", "NS")
    }
  
  # Create rest of column dynamically with comparison table
  for (i in seq_len(nrow(comparisons))) {
    ctrl <- comparisons$ctrl[i]
    case <- comparisons$cases[i]
    
    comp_id <- paste0(case, "_vs_", ctrl)
    
    # Store significants site per test
    sites_welch_i <- comparison_results[[comp_id]][["Welch"]]
    sites_wilcox_i <- comparison_results[[comp_id]][["Wilcoxon"]]
    
    # Add Welch column 
    colname_welch <- paste0(comp_id, "_Welch")
    results[[colname_welch]] <- ifelse(all_sites %in% sites_welch_i, "Significant", "NS")
    
    # Add Wilcoxon column
    colname_wilcox <- paste0(comp_id, "_Wilcoxon")
    results[[colname_wilcox]] <- ifelse(all_sites %in% sites_wilcox_i, "Significant", "NS")
  }
  
  # Transform in data frame
  summary_table <- as.data.frame(results)
}