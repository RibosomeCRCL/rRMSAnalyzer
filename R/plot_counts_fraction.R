#' Plot count distribution among RNAs for each sample
#'
#' @param ribo A RiboClass object.
#'
#' @return a ggplot object
#' @export
#'
#' @examples
#' data("ribo_toy")
#' plot_counts_fraction(ribo_toy)
plot_counts_fraction <- function(ribo) {

  all_data <- dplyr::bind_rows(ribo$data, .id = "samplesid")
  
  all_data_sums <- all_data |>
    dplyr::group_by(samplesid, rna) |>
    dplyr::summarise(sum.counts = sum(count, na.rm = T),
                     n0 = length(which(count < 5)))
  
  all_data_sums$rna <- factor(all_data_sums$rna, levels = ribo$rna_names$current_name)
  
  all_data_sums$samplesid <- factor(all_data_sums$samplesid, levels = names(ribo[["data"]]))
  
  counts_plot <- ggplot2::ggplot(all_data_sums, aes(x = samplesid, y = sum.counts, fill = rna))
  
  counts_plot <- counts_plot + ggplot2::geom_bar(stat = "identity",
                                                 position = "fill") +
    ggplot2::scale_fill_brewer(palette = "Spectral") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = element_text(angle = 45, size = 7, hjust = 1)) +
    labs(title = "Counts distribution per sample", 
         x = "Sample",
         y = "Counts fraction")
  return(counts_plot)
}
