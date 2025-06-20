#' Plot counts for a given position on a given RNA
#' 
#'
#' @param ribo A RiboClass object
#' @param rna Name of RNA where the position is located
#' @param pos Position on RNA on which the view will be centered
#' @param samples Samples to display. "all" will display all samples.
#' @param flanking Number of sites to display on the left/right of the selected position.
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' data("ribo_toy")
#' ribo_toy <- rename_rna(ribo = ribo_toy)
#' plot_counts_env(ribo = ribo_toy, rna = "5.8S", pos = 15)
plot_counts_env <- function(ribo = NULL, rna = NULL, pos = NULL, samples = "all", flanking = 6, condition = NULL) {
  #in this code condition = column name of the metadata BUT [[condition]] and .data[[condition]] = value taken by the parameter "condition" in the function above
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
  count_data <- extract_data(ribo = ribo, col = "count") # to format the data
  
  # extract metadata to merge it with count_data
  metadata <- ribo$metadata[, c("samplename", condition), drop = FALSE]
  #metadata <- ribo$metadata[, c("samplename", "condition", .data[[condition]] ) ] 
  
  # extract the information of position
  which_pos <- which(count_data$named_position == pos_of_interest)
  
  # positions around the position of interest
  which_positions <- c((which_pos - flanking):which_pos, (which_pos + 1):(which_pos + flanking))
  
  #filter to have just the count for the flanking region
  count_data <- count_data[which_positions,] 
  
  # add column for position without RNA before
  count_data$new_position <- c((pos - flanking):pos, (pos + 1):(pos + flanking)) 
  
  # Fusion with metadata to have condition information
  count_transform <- tidyr::gather(count_data[,-1], "samplename", "count", -new_position) 
  
  # Add of new_position and condition column specified
  count_transform <- merge(count_transform, metadata, by = "samplename", all.y = FALSE)  
  
  # check if there are other modifications in the window
  other_mod <- which(!is.na(ribo$data[[1]][which_positions[-ceiling(length(which_positions) / 2)],"site"]))
  
  other_mod_pos <- count_data$new_position[other_mod]
  
  # Count sample number by condition
  count_transform <- count_transform %>%
    dplyr::group_by(if (!is.null(condition)) .data[[condition]] else new_position)
  
  # Count samples number per condition or new_position if there is no condition
  sample_counts <- count_transform %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop")
  
  # Create name condition
  if (!is.null(condition)) {
    labels_with_counts <- setNames(
    paste0(sample_counts[[condition]], " (n=", sample_counts$n, ")"), 
    sample_counts[[condition]]
  )
  }
  # 4 ggplots  
  # -----------------------------------1-------------------------------------------
  if (samples[1] == "all" && is.null(condition)) { #if sample = all and no condition specified  
    
    plot_to_return <- ggplot(data = count_transform) +
      geom_boxplot(aes(x = new_position, y = log10(count), group = new_position)) +
      geom_boxplot(data = count_transform[which(count_transform$new_position %in% other_mod_pos),], 
                   aes(x = new_position, y = log10(count), group = new_position), 
                   fill = "pink",
                   width = 0.8) +
      geom_boxplot(data = count_transform[which(count_transform$new_position == pos),], 
                   aes(x = new_position, y = log10(count), group = new_position), 
                   fill = "green",
                   width = 0.8) +
      theme_bw() +
      labs(title = paste("Count profile for",length(ribo[["data"]]),"samples"),
           subtitle = paste("RNA:",rna),
           y = "log10(count)",
           x = "Position") +
      scale_x_continuous(labels = min(count_transform$new_position):max(count_transform$new_position), 
                         breaks = min(count_transform$new_position):max(count_transform$new_position)) +
      scale_y_continuous(limits = c(0, NA)) +
      #highligth the minimum coverage
      geom_hline(yintercept = 2 , 
                 linewidth = 1,
                 linetype = "11",
                 color = "red") +
      annotate("text",
               label = "Coverage limit",
               x = pos - flanking,
               y = 2 / 1.02,
               color = "red") +
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
  
  # -----------------------------------2-------------------------------------------  
  if ((all(samples %in% names(ribo[["data"]]))) && is.null(condition)) { 
    
    plot_to_return <- ggplot(data = count_transform, aes(x = new_position, y = log10(count), group = samplename)) + 
      geom_point(size = 3) +
      geom_line(aes(col = samplename), linewidth = 2) + 
      theme_bw() + 
      labs(title = paste("Count profile for",length(samples),"samples"),
           subtitle = paste("RNA:",rna),
           y = "log10(count)",
           x = "position") +
      scale_x_continuous(labels = min(count_transform$new_position):max(count_transform$new_position), 
                         breaks = min(count_transform$new_position):max(count_transform$new_position)) +
      scale_y_continuous(limits = c(0, NA)) +
      #highligth the minimum coverage
      geom_hline(yintercept = 2 , 
                 linewidth = 1,
                 linetype = "11",
                 color = "red") +
      annotate("text",
               label = "Coverage limit",
               x = pos - flanking,
               y = 2 / 1.02,
               color = "red") +
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
                 color = "pink") + 
      geom_vline(xintercept = pos,
                 linewidth = 1,
                 linetype = "11",
                 color = "green")
    
  }
  
  # -----------------------------------3-------------------------------------------
  
  if (samples[1] == "all" && !is.null(condition)) { # if samples = all and condition specified

    # Extract uniques modalities of condition variable
    modalities <- unique(count_transform[[condition]])

    #Verifying that there is just two modality
    if (length(modalities) != 2) { #Verify that there is only 2 modalities in the condition
       stop("Variable condition needs to have only 2 modalities")
     }

    # Compute medians per modalities
    median_mod1 <- median(log10(count_transform$count[count_transform[[condition]] == modalities[1]]), na.rm = TRUE)
    median_mod2 <- median(log10(count_transform$count[count_transform[[condition]] == modalities[2]]), na.rm = TRUE)


    plot_to_return <- ggplot(data = count_transform) +
      aes(x = new_position,
          y = log10(count),
          color = .data[[condition]], #for shape color
          fill = NA, #no fill color
          group = interaction(new_position, .data[[condition]])) + #to group by position and specified condition

      geom_boxplot(position = position_dodge(width = 0.7),  #to avoid overlap
                   width = 0.6,
                   fill = NA) +

      # do boxplot filled for specified position and other modifications in the window
      geom_boxplot(data = count_transform[which(count_transform$new_position %in% other_mod_pos),], #specified position
                   aes(x = new_position, y = log10(count), group = interaction(new_position, .data[[condition]])),
                   position = position_dodge(width = 0.7),
                   fill = "pink",
                   alpha = 0.5,
                   width = 0.6,
                   show.legend = FALSE) +

      geom_boxplot(data = count_transform[which(count_transform$new_position == pos),], # other modification
                   aes(x = new_position, y = log10(count), group = interaction(new_position, .data[[condition]])),
                   position = position_dodge(width = 0.7),
                   fill = "green",
                   alpha = 0.5,
                   width = 0.6,
                   show.legend = FALSE) +
      scale_color_manual(values=c("red", "#3182bd"),,
                         labels = labels_with_counts) +

      theme_bw() +

      labs(title = paste("Count profile for", length(ribo[["data"]]), "samples"),
           subtitle = paste("RNA:", rna),
           y = "log10(count)",
           x = "Position") +

      scale_x_continuous(labels = min(count_transform$new_position):max(count_transform$new_position), # to have all the position in x axis
                         breaks = min(count_transform$new_position):max(count_transform$new_position)) +
      scale_y_continuous(limits = c(0, NA)) +
      #highligth the minimum coverage
      geom_hline(yintercept = 2 ,
                 linewidth = 1,
                 linetype = "11",
                 color = "#fc9272") +
      annotate("text",
               label = "Coverage limit",
               x = pos - flanking,
               y = 2 / 1.02,
               color = "#fc9272") +

      # Add madians for each modality of condition
      geom_hline(yintercept = median_mod1, linewidth = 1, linetype = "11", color = "red") +
      annotate("text", label = paste("Counts median", modalities[1]),
               x = pos - flanking + 1, y = median_mod1 / 0.985, color = "red") +

      geom_hline(yintercept = median_mod2, linewidth = 1, linetype = "11", color = "#3182bd") +
      annotate("text", label = paste("Counts median", modalities[2]),
               x = pos - flanking + 1, y = median_mod2 / 0.985, color = "#3182bd")
  }

  # -----------------------------------4-------------------------------------------
  
  if ((all(samples %in% names(ribo[["data"]]))) && !is.null(condition)) { 
    # Extract uniques modalities of condition variable
    modalities <- unique(count_transform$condition)
    
    # Verifying that there is just two modality
    if (length(modalities) != 2) { 
      stop("Variable condition needs to have only 2 modalities")
    }
    
    # Compute madians per modalities
    median_mod1 <- median(log10(count_transform$count[count_transform[[condition]] == modalities[1]]), na.rm = TRUE)
    median_mod2 <- median(log10(count_transform$count[count_transform[[condition]] == modalities[2]]), na.rm = TRUE)
  
    plot_to_return <- ggplot(data = count_transform) + 
      aes(x = new_position,  
          y = log10(count), 
          color = .data[[condition]], #for shape color
          fill = NA, #no fill color
          group = interaction(new_position, .data[[condition]])) + #to group by position and specified condition
      
      geom_boxplot(position = position_dodge(width = 0.7),  #to avoid overlap 
                   width = 0.6,
                   fill = NA) +  
      
      # do boxplot filled for specified position and other modifications in the window
      geom_boxplot(data = count_transform[which(count_transform$new_position %in% other_mod_pos),], #specified position
                   aes(x = new_position, y = log10(count), group = interaction(new_position, .data[[condition]])), 
                   position = position_dodge(width = 0.7),
                   fill = "pink",
                   alpha = 0.5,
                   width = 0.6,
                   show.legend = FALSE) +
      
      geom_boxplot(data = count_transform[which(count_transform$new_position == pos),], # other modification
                   aes(x = new_position, y = log10(count), group = interaction(new_position, .data[[condition]])), 
                   position = position_dodge(width = 0.7),
                   fill = "green",
                   alpha = 0.5,
                   width = 0.6,
                   show.legend = FALSE) +
      
      theme_bw() +
      
      labs(title = paste("Count profile for", length(ribo[["data"]]), "samples"),
           subtitle = paste("RNA:", rna),
           y = "log10(count)",
           x = "Position") +
      
      scale_x_continuous(labels = min(count_transform$new_position):max(count_transform$new_position), # to have all the position in x axis
                         breaks = min(count_transform$new_position):max(count_transform$new_position)) + 
      scale_y_continuous(limits = c(0, NA)) +

      #highligth the minimum coverage
      geom_hline(yintercept = 2 , 
                 linewidth = 1,
                 linetype = "11",
                 color = "#fc9272") +
      annotate("text",
               label = "Coverage limit",
               x = pos - flanking,
               y = 2 / 1.02,
               color = "#fc9272") +
      
      # Add madians for each modality of condition
      geom_hline(yintercept = median_mod1, linewidth = 1, linetype = "11", color = "red") +  
      annotate("text", label = paste("Counts median", modalities[1]), 
               x = pos - flanking + 1, y = median_mod1 / 0.985, color = "red") + 
      
      geom_hline(yintercept = median_mod2, linewidth = 1, linetype = "11", color = "#3182bd") +  
      annotate("text", label = paste("Counts median", modalities[2]), 
               x = pos - flanking + 1, y = median_mod2 / 0.985, color = "#3182bd") + 
      
      scale_color_manual(values = c("red", "#3182bd"),
                         ,
                         labels = labels_with_counts)
  }

  return(plot_to_return)
}