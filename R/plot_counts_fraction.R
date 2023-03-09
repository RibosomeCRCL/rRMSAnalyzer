
plot_counts_fraction <- function(ribo = NULL) {
  #output
  # faire les check
  all_data <- dplyr::bind_rows(ribo$data, .id = "samplesid")
  
  all_data_sums <- all_data %>%
    group_by(samplesid, rna) %>%
    summarise(sum.counts = sum(count, na.rm = T), n0 = length(which(count < 5)))
  
  all_data_sums$rna <- factor(all_data_sums$rna, levels = ribo$rna_names$current_name)
  
  all_data_sums$samplesid <- factor(all_data_sums$samplesid, levels = names(ribo[["data"]]))
  
  counts_plot <- ggplot(all_data_sums, aes(x = samplesid, y = sum.counts, fill = rna))
  
  counts_plot <- counts_plot + geom_bar(stat = "identity", position = "fill") +
    scale_fill_brewer(palette = "Spectral") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, size = 7, hjust = 1)) +
    labs(title = "Counts distribution per sample", 
         x = "Sample",
         y = "Counts fraction")
  return(counts_plot)
}
