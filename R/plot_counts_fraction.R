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
  # NSE fix
  samplesid <- rna <- count <- sum.counts <- NULL

  all_data <- dplyr::bind_rows(ribo$data, .id = "samplesid")

  all_data_sums <- all_data |>
    dplyr::group_by(samplesid, rna) |>
    dplyr::summarise(sum.counts = sum(count, na.rm = TRUE),
                     n0 = length(which(count < 5)))

  all_data_sums$rna <- factor(all_data_sums$rna, levels = ribo$rna_names$current_name)

  all_data_sums$samplesid <- factor(all_data_sums$samplesid, levels = names(ribo[["data"]]))

  counts_plot <- ggplot2::ggplot(all_data_sums, aes(x = samplesid, y = sum.counts, fill = rna))

  counts_plot <- counts_plot + ggplot2::geom_bar(stat = "identity",
                                                 position = "fill") +
    ggplot2::scale_fill_brewer(palette = "Spectral") +
    ggplot2::theme_bw() +
    #add théorical mean fraction 
    ggplot2::geom_hline(yintercept = 0.60381342, linewidth = 1, linetype = "11", color = "blue") +  
    ggplot2::annotate("text", label = "28S",
                      levels(factor(all_data_sums$samplesid))[length(levels(factor(all_data_sums$samplesid)))],
                      y = 0.604 + 0.02, color = "blue") + 
    
    ggplot2::geom_hline(yintercept = 0.9250826, linewidth = 1, linetype = "11", color = "darkgreen") +  #0.604+0.32126916
    ggplot2::annotate("text", label = "18S",
                      levels(factor(all_data_sums$samplesid))[length(levels(factor(all_data_sums$samplesid)))],
                      y = 0.925 + 0.02, color = "darkgreen") + 
    
    ggplot2::geom_hline(yintercept = 0.9711803, linewidth = 1, linetype = "11", color = "orange2") +  #0.9250826+0.04609767
    ggplot2::annotate("text", label = "5.8S",
                      x = levels(factor(all_data_sums$samplesid))[length(levels(factor(all_data_sums$samplesid)))], 
                      y = 0.9711 + 0.02, color = "black") + 
    
    ggplot2::geom_hline(yintercept = 1, linewidth = 1, linetype = "11", color = "red") +  #0.9711803+ 0.02881975
    ggplot2::annotate("text", label = "5S",
                      x = levels(factor(all_data_sums$samplesid))[1], 
                      y = 0.9995 , color = "black") + 
    
    ggplot2::theme(axis.text.x = element_text(angle = 90, size = 7, hjust = 1)) +
    labs(title = "Counts distribution per sample",
         x = "Sample",
         y = "Counts fraction")
  return(counts_plot)
}