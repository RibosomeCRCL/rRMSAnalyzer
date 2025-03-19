#' Plot statitic
#'
#' @param ribo a RiboClass
#' @param site sites to plot
#' @param res_pv df of p values extracted from res_pv.R
#' @param pthr p value threshold
#'
#' @returns a boxplot
#' @export
#'
#' @examples plot_stat(data = ribo, site = NULL, res_pv = res_pv, pthr = 0.05)
plot_stat <- function(ribo = ribo_adj_annot, site = NULL, res_pv = res_pv, pthr = NULL) {
  library(ggplot2)
  
    # Cas où on ne peut pas afficher de plot
  if (is.null(site) && (is.null(res_pv) || is.null(pthr))) {
    stop("Erreur : Vous devez spécifier au moins un site")
  }
  if (!is.null(site) && is.null(res_pv) && !is.null(pthr)) {
    stop("Erreur : res_pv est requis si pthr est fourni")
  }

  # filter sites
  if (!is.null(res_pv) && !is.null(pthr)) {
    significant_sites <- res_pv %>% dplyr::filter(p_adj < pthr) %>% dplyr::pull(annotated_sites)
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
  stat <- ggplot(data, aes(x = .data[[condition_col]], y = c_score, fill = .data[[condition_col]])) +
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
