#' Plot comparison of mediane of c-score between two conditions
#'
#' @param ribo A RiboClass object.
#' @param condition_col condition column for the plot
#' @param ech modality 1 of condition_col
#' @param ctrl modality 2 of condition_col
#' 
#' @returns two ggplot objects
#' @export
#'
#' @examples
#' plot_comparison_medianes(ribo_adj_annot, params$condition_col, params$ech, params$ctrl)

plot_comparison_medianes <- function(ribo = NULL, condition_col = NULL,  ech = NULL, ctrl = NULL) {
  
  ribo_df <- extract_data(ribo_adj_annot,
                        col = "cscore",
                        only_annotated = TRUE)
# 7.1. transforming df

cond1 <- ribo$metadata %>% #filter sample of your 1st condition 
  dplyr::filter(!!sym(condition_col) == paste0(ech)) 

cond2 <- ribo$metadata %>% #filter sample of your 2nd condition 
  dplyr::filter(!!sym(condition_col) == paste0(ctrl)) 

cond1_2 <- unique(c(cond1$samplename, cond2$samplename)) # merge both of your condition into a variable 

new_metadata_sites <- ribo_df %>%
  tidyr::pivot_longer(cols = -site,         # All columns but "site"
               names_to = "samplename", 
               values_to = "cscore") %>%
  dplyr::filter(samplename %in% cond1_2)  # Filter only sample present in cond1 et cond2


metadata_data_total <- new_metadata_sites %>%
  dplyr::group_by(site) %>%
  dplyr::summarize(median_Cscore_metadata_ech_ctrl = median(cscore[samplename %in% cond1$samplename]), # take samplename into cond1
            median_Cscore_metadata_ech = median(cscore[samplename %in% cond2$samplename]), # take samplename into cond2
            difference = median_Cscore_metadata_ech - median_Cscore_metadata_ech_ctrl) # compute the difference between both



metadata_data_total <- metadata_data_total %>%
  dplyr::mutate(site = as.character(site))

# keep site order (to debug)
# metadata_data_total <- metadata_data_total %>%
#   dplyr::mutate(site = factor(site, levels = rev(unique(human_methylated$Nomenclature))))

#create difference_cat column
metadata_data_total$difference_cat <- ifelse(abs(metadata_data_total$difference) < 0.05, "No difference",
                                             ifelse(metadata_data_total$difference > 0.05, "Increase",
                                                    ifelse(metadata_data_total$difference < 0.05, "Decrease", "ok")))

# metadata_data_total <- metadata_data_total %>%
#   dplyr::mutate(
#     difference_cat = dplyr::case_when(
#       abs(difference) < 0.05 ~ "No difference",
#       difference >= 0.05     ~ "Increase",
#       difference <= -0.05    ~ "Decrease"
#     ),
#     difference_cat = factor(difference_cat, levels = c("Increase", "No difference", "Decrease"))
#   )

# create a color vector for tags according to difference_cat
site_colors <- ifelse(metadata_data_total$difference_cat == "Increase", "red",
                      ifelse(metadata_data_total$difference_cat == "Decrease", "blue", "black"))
# # /!\ Vectorized input to `element_text()` is not officially supported.
# ℹ Results may be unexpected or may change in future versions of ggplot2. 

#site_colors <- c("Increase" = "red", "No difference" = "black", "Decrease" = "blue")

#7.3. Graphique n°1: vizualisation of absolute cscore value
part1 <- ggplot(metadata_data_total, aes(x=site, y=difference, label=site)) +
  geom_bar(stat='identity', aes(fill=difference_cat), width = .9) +
  scale_fill_manual(name="Trend of 2'Ome variation",
                    labels = c("Decrease", "Increase", "No difference"),
                    values = c("Increase"="red", "Decrease"="blue", "No difference" = "grey"))  +
  labs(title= "Representation of differential median C-score") +
  geom_hline(yintercept = 0.05, linetype = "11") +
  geom_hline(yintercept = -0.05, linetype = "11") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8, color = site_colors), # tag size and color of y axis
        legend.position = "bottom",
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.minor = element_blank()) +
  annotate("rect", ymin = -0.05, ymax = 0.05,
           xmin = 0, xmax = 113, fill = "darkgrey", alpha = 0.2) +
  xlab("2'Ome sites") +
  ylab(bquote("\u0394"~ "median C-score" ~ "(" * .(ech) ~ "-" ~ .(ctrl) * ")")) +
  coord_flip() +
  scale_y_continuous(sec.axis = sec_axis(trans = ~ ., name = NULL, breaks = NULL), limits = c(-1, 1))

#7.4. Graphique n°2: difference of cscore vizualisation 
part2 <- ggplot(metadata_data_total) +
  geom_point(data = metadata_data_total[,c(1,3)],
             aes(x = site, y = median_Cscore_metadata_ech, color = paste0(ech)), 
             size = 2) +
  geom_point(data = metadata_data_total[,c(1,2)],
             aes(x = site, y = median_Cscore_metadata_ech_ctrl, color = paste0(ctrl)),
             size = 2) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.minor = element_blank()) +
  geom_segment(data = metadata_data_total[which(metadata_data_total$difference_cat %in% c("Decrease", "Increase")), ],
               aes(x = site, y = median_Cscore_metadata_ech, xend = site, yend = median_Cscore_metadata_ech_ctrl, color = difference_cat),
               alpha = 1,
               linewidth = 1) +
  coord_flip() +
  scale_y_continuous(name = "Median c-score",sec.axis = sec_axis(trans = ~ ., name = NULL, breaks = NULL)) +
  scale_color_manual(
    values = setNames(c("orchid", "seagreen", "grey", "blue", "red"),
                      c(ctrl, ech, "No difference", "Decrease", "Increase")), # format colours dynamically 
    name = "Condition",) +
  theme(legend.position = "bottom")
  
  # 7.5. combine graphs
  plot_comparison_mediane <- patchwork::wrap_plots(part1, part2, widths = c(2,1))

  return(plot_comparison_mediane)
}
