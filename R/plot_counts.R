#' Plot counts for a given position on a given RNA
#' 
#'
#' @param ribo A RiboClass object
#' @param rna Name of RNA where the position is located
#' @param pos Position on RNA on which the view will be centered
#' @param samples Samples to display. "all" will display all samples.
#' @param flanking Number of sites to display on the left/right of the selected position.
#'
#' @return
#' @export
#'
#' @examples
#' data("ribo_toy")
#' ribo_toy <- rename_rna(ribo = ribo_toy)
#' plot_counts_env(ribo = ribo_toy, rna = "5.8S", pos = 15)
plot_counts_env <- function(ribo = NULL, rna = NULL, pos = NULL, samples = "all", flanking = 6) {
  new_position <- count <- NULL
  
  #Check for ribo
  if (is.null(ribo)) {stop("MISSING parameter: please provide a RiboClass!")}
  if (!inherits(ribo, "RiboClass")) {stop("ribo argument is not a RiboClass!")}
  
  # check if rna is ok
  if(!(rna %in% ribo[["rna_names"]][["current_name"]])) {
    stop("the RNA names given do not exist in the RiboClass")
  }
  # check for position
  if(!(is.numeric(pos))) {stop("pos should be a number")}
  if (table(ribo$data[[1]]$rna)[rna] < pos) {stop(paste(pos, " is higher than the lenth of", rna))}
  
  # check for sample
  
  if(samples[1] == "all") {
    ribo <- ribo
  } else if (all(samples %in% names(ribo[["data"]]))) {
    ribo <- keep_ribo_samples(ribo = ribo, samples_to_keep = samples)
  } else {
    stop("Samples name should be from ", toString(names(ribo[["data"]])))
  }
  
  # check for flanking
  if(!(is.numeric(flanking))) {stop("flanking should be a number")}
  
  
  # format position of interest according to named_position of RiboClass
  pos_of_interest <- paste(rna, formatC(pos,width = 4,flag = "0"), sep = "_")
  
  # Extract count data from RiboClass
  count_data <- extract_data(ribo = ribo, col = "count")
  
  # extract the information
  which_pos <- which(count_data$named_position == pos_of_interest)
  # positions around the position of interest
  which_positions <- c((which_pos - flanking):which_pos, (which_pos + 1):(which_pos + flanking))
  
  count_data <- count_data[which_positions,]
  
  count_data$new_position <- c((pos - flanking):pos, (pos + 1):(pos + flanking))
  
  count_transform <- tidyr::gather(count_data[,-1], "samples", "count", -new_position)
  
  
  # check if there are other modifications in the window
  
  other_mod <- which(!is.na(ribo$data[[1]][which_positions[-ceiling(length(which_positions) / 2)],"site"]))
  
  other_mod_pos <- count_data$new_position[other_mod]
  
  
  # ggplot 
  
  if (samples[1] == "all") {
    
    plot_to_return <- ggplot(data = count_transform) +
      geom_boxplot(aes(x = new_position, y = log10(count), group = new_position)) +
      geom_boxplot(data = count_transform[which(count_transform$new_position %in% other_mod_pos),], 
                   aes(x = new_position, y = log10(count), group = new_position), 
                   fill = "lightblue",
                   width = 0.8) +
      geom_boxplot(data = count_transform[which(count_transform$new_position == pos),], 
                   aes(x = new_position, y = log10(count), group = new_position), 
                   fill = "lightgreen",
                   width = 0.8) +
      theme_bw() +
      labs(title = paste("Count profile for",length(ribo[["data"]]),"samples"),
           subtitle = paste("RNA:",rna),
           y = "log10(count)",
           x = "Position") +
      scale_x_continuous(labels = min(count_transform$new_position):max(count_transform$new_position), 
                         breaks = min(count_transform$new_position):max(count_transform$new_position)) +
      {if (min(count_transform$count) < 100) geom_hline(yintercept = 2 ,
                                                        linewidth = 1,
                                                        linetype = "11",
                                                        color = "red")} +
      {if (min(count_transform$count) < 100) annotate("text", 
                                                      label = "Coverage limit",
                                                      x = pos - flanking,
                                                      y = 2 / 1.02,
                                                      color = "red")} +
      geom_hline(yintercept = stats::median(log10(count_transform$count) , na.rm = TRUE),
                 linewidth = 1,
                 linetype = "11",
                 color = "darkred") +
      annotate("text", 
               label = "Counts median",
               x = pos - flanking,
               y = stats::median(log10(count_transform$count) / 0.985, na.rm = TRUE),
               color = "darkred")
      

  } 
  if ((all(samples %in% names(ribo[["data"]])))) {
    
    plot_to_return <- ggplot(data = count_transform, aes(x = new_position, y = log10(count), group = samples)) +
      geom_point(size = 3) +
      geom_line(aes(col = samples), linewidth = 2) +
      theme_bw() + 
      labs(title = paste("Count profile for",length(samples),"samples"),
           subtitle = paste("RNA:",rna),
           y = "log10(count)",
           x = "position") +
      scale_x_continuous(labels = min(count_transform$new_position):max(count_transform$new_position), 
                         breaks = min(count_transform$new_position):max(count_transform$new_position)) +
      {if (min(count_transform$count) < 100) geom_hline(yintercept = 2 ,
                                                        linewidth = 1,
                                                        linetype = "11",
                                                        color = "red")} +
      {if (min(count_transform$count) < 100) annotate("text", 
                                                      label = "Coverage limit",
                                                      x = pos - flanking + 1,
                                                      y = 2 / 1.02,
                                                      color = "red")} +
      geom_hline(yintercept = stats::median(log10(count_transform$count) , na.rm = TRUE),
                 linewidth = 1,
                 linetype = "11",
                 color = "darkred") +
      annotate("text", 
               label = "Counts median",
               x = pos - flanking + 1,
               y = stats::median(log10(count_transform$count) / 0.985, na.rm = TRUE),
               color = "darkred") + 
      geom_vline(xintercept = other_mod_pos,
                 linewidth = 1,
                 linetype = "11",
                 color = "lightblue") + 
      geom_vline(xintercept = pos,
                 linewidth = 1,
                 linetype = "11",
                 color = "lightgreen")
      
  }
  
  return(plot_to_return)

}


# rajoute une ligne à 100
