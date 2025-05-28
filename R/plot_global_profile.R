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

plot_global_profile <- function(ribo = ribo_adj_annot, condition_col = NULL) { #, ech = NULL, ctrl = NULL) {
  library(rlang)
  #Check for ribo
  if (is.null(ribo)) {stop("MISSING parameter: please provide a RiboClass!")}
  if (!inherits(ribo, "RiboClass")) {stop("ribo argument is not a RiboClass!")}
  
  metadata_x <- ribo$metadata # extract metadata from ribo object
  ribom <- extract_data(ribo, col = "cscore", only_annotated = TRUE, position_to_rownames = TRUE)
  
  # Transform ribom into long format excluding -annotated_sites
  ribom_long <- ribom %>% 
    tibble::rownames_to_column("annotated_sites") %>%  # Add annotated_sites as column
    tidyr::pivot_longer(cols = -annotated_sites, names_to = "samplename", values_to = "c_score")  # keep annotated_sites outside pivoting

# fusion with metadata
  ribom_merged <- ribom_long %>% 
    dplyr::left_join(metadata_x, by = "samplename") # left join between metadata and ribom_long 
  
  result <- ribom_merged %>%
    dplyr::group_by(annotated_sites, !!sym(condition_col)) %>%
    dplyr::summarise(mean_c_score = mean(c_score, na.rm = TRUE),
              sd_c_score = sd(c_score, na.rm = TRUE),
              .groups = "drop") #%>%
    #tidyr::pivot_wider(names_from = !!sym(condition_col), values_from = c(mean_c_score, sd_c_score))
  
  #keep order as in human methylated 
  result <- result %>% 
    dplyr::mutate(annotated_sites = factor(annotated_sites, 
                                    levels = unique(human_methylated$Nomenclature))) %>%
    dplyr::arrange(annotated_sites)

  # -----------------------------------1----------------------------------------
  #if (length(unique(condition_col)) <= 4) {   
    # drawing the graph
#     global_profile <- ggplot(result, aes(x = annotated_sites, group = 1)) +
#       
#       geom_line(aes(y = .data[[paste0("mean_c_score_", ech)]], color = paste0(ech))) + 
#       geom_line(aes(y = .data[[paste0("mean_c_score_", ctrl)]], color = paste0(ctrl))) + 
#       
#       geom_ribbon(aes(ymin = .data[[paste0("mean_c_score_", ctrl)]] - .data[[paste0("sd_c_score_", ech)]], 
#                       ymax = .data[[paste0("mean_c_score_", ctrl)]] + .data[[paste0("sd_c_score_", ech)]], 
#                       fill = paste0(ech)), alpha = 0.3) +
#       
#       geom_ribbon(aes(ymin = .data[[paste0("mean_c_score_", ctrl)]] - .data[[paste0("sd_c_score_", ctrl)]], 
#                       ymax = .data[[paste0("mean_c_score_", ctrl)]] + .data[[paste0("sd_c_score_", ctrl)]],
#                       fill = paste0(ctrl)), alpha = 0.3) +
#       
#       scale_x_discrete(guide = guide_axis(angle = 90)) +
#       scale_y_continuous(limits = c(0, 1)) +
#       labs(x = "RNA 2'Ome Methylation sites", y = "C-score", 
#            color = "Condition", fill = "Condition",
#            title = "2'Ome Globale Profile") +
#       theme_bw() +
#       scale_color_brewer(palette = "Paired")   
#     
#   return(global_profile)
#   } else {
#     stop("There is to much condition into the condition to draw a clear plot")
#   }
# }

  # -----------------------------------2----------------------------------------
    # Extraire les colonnes dynamiques
    # mean_cols <- grep("^mean_c_score_", names(result), value = TRUE)
    # sd_cols   <- grep("^sd_c_score_", names(result), value = TRUE)
    # conditions <- sub("mean_c_score_", "", mean_cols)
    
  if (length(unique(condition_col)) > 4) stop("Too many conditions to draw a clear plot.")
    
    global_profile <- ggplot(result, aes(x = annotated_sites, y = mean_c_score , group = condition , color = condition, fill = condition)) +
      geom_line(linewidth = 1) +
      geom_ribbon(aes(ymin = mean_c_score - sd_c_score, ymax = mean_c_score + sd_c_score), alpha = 0.3) +
      scale_x_discrete(guide = guide_axis(angle = 90)) +
      scale_y_continuous(limits = c(0, 1)) +
      labs(
        x = "RNA 2'Ome Methylation sites", y = "C-score",
        color = "Condition", fill = "Condition",
        title = "2'Ome Global Profile"
      ) +
      theme_bw() +
      scale_color_brewer(palette = "Accent") +
      scale_fill_brewer(palette = "Accent")
    
    return(global_profile)
  }
