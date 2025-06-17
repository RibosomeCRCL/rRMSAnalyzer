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
  dplyr::summarize(median_Cscore_metadata_ech = median(cscore[samplename %in% cond1$samplename]), # take samplename into cond1
            median_Cscore_metadata_ech_ctrl = median(cscore[samplename %in% cond2$samplename]), # take samplename into cond2
            difference = median_Cscore_metadata_ech - median_Cscore_metadata_ech_ctrl) # compute the difference between both



metadata_data_total <- metadata_data_total %>%
  dplyr::mutate(site = as.character(site))

# keep site order (doesn't work yet)
# metadata_data_total <- metadata_data_total %>%
#   dplyr::mutate(site = factor(site, levels = rev(unique(human_methylated$Nomenclature))))

#create difference_cat column
metadata_data_total$difference_cat <- ifelse(abs(metadata_data_total$difference) < 0.05, "No difference",
                                             ifelse(metadata_data_total$difference > 0.05, "Increase",
                                                    ifelse(metadata_data_total$difference < 0.05, "Decrease", "ok")))

# Force level of color to keep legend every time 

metadata_data_total$difference_cat <- factor(metadata_data_total$difference_cat,
                                             levels = c("Decrease", "Increase", "No difference"))

# create a color vector for tags according to difference_cat
site_colors <- ifelse(metadata_data_total$difference_cat == "Increase", "red",
                      ifelse(metadata_data_total$difference_cat == "Decrease", "blue", "black"))

#site_colors <- c("Increase" = "red", "No difference" = "black", "Decrease" = "blue")

# add column for colouring x axis 
metadata_data_total$site_color <- site_colors  

# Force order level to don't disrupt legend
metadata_data_total$legend_color <- factor(metadata_data_total$difference_cat,
                                           levels = c("Decrease", "Increase", "No difference"))

#7.3. Graphique n°1: vizualisation of absolute cscore value
part1 <- ggplot(metadata_data_total, aes(x = site, y = difference, label = site)) +
  geom_bar(stat = 'identity', aes(fill = difference_cat), width = .9) +
  scale_fill_manual(
    name = "Trend of 2'Ome variation",
    labels = c("Decrease", "Increase", "No difference"),
    values = c("Increase" = "red", "Decrease" = "blue", "No difference" = "grey"),
    drop = FALSE
  ) +
  geom_hline(yintercept = 0.05, linetype = "11") +
  geom_hline(yintercept = -0.05, linetype = "11") +
  annotate("rect", ymin = -0.05, ymax = 0.05, xmin = 0, xmax = 113, fill = "darkgrey", alpha = 0.2) +
  labs(title = "Representation of differential median C-score") +
  theme_bw() +
  theme(
    axis.text.y = element_blank(),  # on supprime les labels automatiques
    axis.ticks.y = element_blank(), # on supprime aussi les ticks
    legend.position = "bottom",
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.minor = element_blank()
  ) +
  geom_text(aes(x = site, y = -1.05, label = site, color = site_color), hjust = 1, size = 2.5) +  # add label manually
  scale_color_identity() +  
  xlab("2'Ome sites") +
  ylab(bquote("\u0394"~ "median C-score" ~ "(" * .(ech) ~ "-" ~ .(ctrl) * ")")) +
  coord_flip() +
  scale_y_continuous(
    sec.axis = sec_axis(transform = ~ ., name = NULL, breaks = NULL), limits = c(-1.1, 1))
# part1 <- ggplot(metadata_data_total, aes(x=site, y=difference, label=site)) +
#   geom_bar(stat='identity', aes(fill=difference_cat), width = .9) +
#   scale_fill_manual(name="Trend of 2'Ome variation",
#                     labels = c("Decrease", "Increase", "No difference"),
#                     values = c("Increase"="red", "Decrease"="blue", "No difference" = "grey"))  +
#   labs(title= "Representation of differential median C-score") +
#   geom_hline(yintercept = 0.05, linetype = "11") +
#   geom_hline(yintercept = -0.05, linetype = "11") +
#   theme_bw() +
#   theme(axis.text.y = element_text(size = 8, color = site_colors), # tag size and color of y axis
#         legend.position = "bottom",
#         plot.subtitle = element_text(hjust = 0.5),
#         panel.grid.minor = element_blank()) +
#   annotate("rect", ymin = -0.05, ymax = 0.05,
#            xmin = 0, xmax = 113, fill = "darkgrey", alpha = 0.2) +
#   xlab("2'Ome sites") +
#   ylab(bquote("\u0394"~ "median C-score" ~ "(" * .(ech) ~ "-" ~ .(ctrl) * ")")) +
#   coord_flip() +
#   scale_y_continuous(sec.axis = sec_axis(transform = ~ ., name = NULL, breaks = NULL), limits = c(-1, 1))

#create a column to tag ech and column (avoid confusion in legend color)
## Store name dynamically 
# condition_labels <- c("ctrl" = ctrl, "ech" = ech)
# 
# df_ech <- metadata_data_total[, c("site", "median_Cscore_metadata_ech")]
# df_ech$Condition <- "ech"
# 
# df_ctrl <- metadata_data_total[, c("site", "median_Cscore_metadata_ech_ctrl")]
# df_ctrl$Condition <- "ctrl"
# 
# # Rename and fusion
# colnames(df_ech)[2] <- "median_cscore"
# colnames(df_ctrl)[2] <- "median_cscore"
# df_all <- rbind(df_ech, df_ctrl)

#7.4. Graphique n°2: difference of cscore vizualisation 
part2 <- ggplot(metadata_data_total) +
  geom_point(data = metadata_data_total[,c(1,2)],
             aes(x = site, y = median_Cscore_metadata_ech, color = paste0(ech)),
             size = 2) +
  geom_point(data = metadata_data_total[,c(1,3)],
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
  scale_y_continuous(name = "Median c-score",sec.axis = sec_axis(transform = ~ ., name = NULL, breaks = NULL)) +
  scale_color_manual(
    values = setNames(c("black", "orange", "grey", "blue", "red"),
                      c(ctrl, ech, "No difference", "Decrease", "Increase")), # format colours dynamically
    name = "Condition",
    drop = FALSE) +
  theme(legend.position = "bottom")
  



# part2 <- ggplot() +
#   geom_point(data = df_all,
#              aes(x = site, y = median_cscore, color = Condition),
#              size = 2) +
#   geom_segment(data = metadata_data_total[metadata_data_total$difference_cat %in% c("Decrease", "Increase"), ],
#                aes(x = site, y = median_Cscore_metadata_ech, xend = site, yend = median_Cscore_metadata_ech_ctrl,
#                    color = difference_cat),
#                alpha = 1,
#                linewidth = 1) +
#   scale_color_manual(
#     values = c("ech" = "orchid", "ctrl" = "seagreen",
#                "Increase" = "red", "Decrease" = "blue", "No difference" = "grey"),
#     labels = c(condition_labels, "Decrease" = "Decrease", "Increase" = "Increase", "No difference" = "No difference"),
#     breaks = c("ctrl", "ech", "Decrease", "Increase", "No difference"),
#     name = "Condition"
#   ) +
#   theme_bw() +
#   theme(axis.title.y = element_blank(),
#         axis.text.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         plot.subtitle = element_text(hjust = 0.5),
#         panel.grid.minor = element_blank(),
#         legend.position = "bottom") +
#   coord_flip() +
#   scale_y_continuous(name = "Median c-score",
#                      sec.axis = sec_axis(transform = ~ ., name = NULL, breaks = NULL))
  # 7.5. combine graphs
  plot_comparison_mediane <- patchwork::wrap_plots(part1, part2, widths = c(2,1))

  return(plot_comparison_mediane)
}
