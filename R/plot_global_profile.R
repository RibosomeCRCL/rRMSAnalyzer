#' Plot global profile for a given condition on a given metadata
#'
#' @param ribo A RiboClass object.
#' @param condition_col condition column for the plot
#' 
#' 
#' @return a ggplot object
#' @export
#' 
#' @examples
#' plot_global_profile(ribo, condition_col, ech, ctrl)

plot_global_profile <- function(ribo = ribo_adj_annot, condition_col = NULL) { 

  #Check for ribo
  if (is.null(ribo)) {stop("MISSING parameter: please provide a RiboClass!")}
  if (!inherits(ribo, "RiboClass")) {stop("ribo argument is not a RiboClass!")}
  
  meta <- ribo$metadata # extract metadata from ribo object
  ribom <- extract_data(ribo, col = "cscore", only_annotated = TRUE, position_to_rownames = TRUE)
  
  # Transform ribom into long format excluding -annotated_sites
  ribom_long <- ribom %>% 
    tibble::rownames_to_column("annotated_sites") %>%  # Add annotated_sites as column
    tidyr::pivot_longer(cols = -annotated_sites, names_to = "samplename", values_to = "c_score")  # keep annotated_sites outside pivoting

# fusion with metadata
  ribom_merged <- ribom_long %>% 
    dplyr::left_join(meta, by = "samplename") # left join between metadata and ribom_long 
  
  result <- ribom_merged %>%
    dplyr::group_by(annotated_sites, !!sym(condition_col)) %>%
    dplyr::summarise(mean_c_score = mean(c_score, na.rm = TRUE),
              sd_c_score = sd(c_score, na.rm = TRUE),
              .groups = "drop") 
  
  #keep order as in human methylated 
  #first order
  fst_order <- (unique(human_methylated$Nomenclature))
  
  # second order
  other_sites <- setdiff(unique(result$annotated_sites), fst_order) # sites only present in result 
  
  #concatenation of both order
  final_order <- c(other_sites, fst_order) #order 5.8S before 18 and 28-S
  
  #use final_order
  result <- result %>%
    dplyr::mutate(annotated_sites = factor(annotated_sites, levels = final_order)) 

  # draw graph
  if (length(unique(condition_col)) > 4) stop("Too many conditions to draw a clear plot.")
    
    global_profile <- ggplot(result, aes(x = annotated_sites, y = mean_c_score , group = !!sym(condition_col) , color = !!sym(condition_col), fill = !!sym(condition_col))) +
      geom_line(linewidth = 1) +
      geom_ribbon(aes(ymin = mean_c_score - sd_c_score, ymax = mean_c_score + sd_c_score), alpha = 0.3) +
      scale_x_discrete(guide = guide_axis(angle = 90)) +
      scale_y_continuous(limits = c(0, 1)) +
      labs(
        x = "rRNA 2'Ome sites", y = "C-score",
        color = condition_col, fill = condition_col,
        title = "Variation in 2'Ome levels between conditions"
      ) +
      theme_bw() +
      scale_color_brewer(palette = "Accent", name = "Condition") +
      scale_fill_brewer(palette = "Accent", name = "Condition")
    
    return(global_profile)
  }
