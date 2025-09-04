#' Plot a boxplot of a RiboClass object counts 
#' @description This plot is useful to check if the samples are alike in their
#' raw counts.
#' @param ribo A RiboClass object.
#' @param color_col Name of the column in the metadata used for coloring samples.
#' @param outlier Show boxplot outlier values.
#' @param horizontal Show boxplot horizontally.
#' @return A ggplot object.
#' @export
#'
#' @examples 
#' data("ribo_toy")
#' boxplot_count(ribo_toy,"run")
 
boxplot_count <- function(ribo, color_col = NA,
                          outlier = TRUE, horizontal = FALSE) {
  ribo_matrix <- extract_data(ribo, "count", position_to_rownames = TRUE)
  
  return(.plot_boxplot_samples(ribo_matrix, "count",
                               ribo[["metadata"]], color_col, outlier,
                               horizontal = horizontal))
}

#' Plot boxplot representing the C-score values of all samples for each
#' individual annotated site
#' Sites are sorted by their median.
#' @inheritParams boxplot_count 
#' @param sort_by Sort sites by median ("median", default) by variance ("var")
#'  or IQR ("iqr").
#' @return a ggplot geom_boxplot
#' @export
#'
#' @examples
#' data("ribo_toy")
#' data("human_methylated")
#' ribo_toy <- rename_rna(ribo_toy)
#' ribo_toy <- annotate_site(ribo_toy,human_methylated)
#' boxplot_cscores(ribo_toy)
#' 
boxplot_cscores <- function(ribo,outlier = TRUE, sort_by = c("median","iqr","var")[1], horizontal = FALSE) { # main function which call inner function .plot_boxplot_sites 
  ribo_m <- extract_data(ribo,only_annotated = TRUE) # extract annotated data
  
  if(nrow(ribo_m) == 0) { # verification that the riboclass is annotated
    stop("No annotated site found. Please use annotate_site() on your RiboClass before calling this function.")
    }
  
  return(.plot_boxplot_sites(ribo_m, # inner function
                      values_to_plot = "cscore",outlier = TRUE, 
                      sort_by = sort_by, # this variable allows to sort the sites annotated according to the median (default) or the IQR or the var
                      horizontal = horizontal))
}

#' Internal function for boxplot_cscore
#'
#' @param matrix Sites x Samples C-score/count matrix (output of extract_data()).
#' @param values_to_plot Value to display in plot.
#' @param outlier Show boxplot outlier values.
#' @param horizontal Show boxplot horizontally.
#'
#' @return A ggplot object.
#' @keywords internal
#' These two functions above and below allow to plot boxplots to visualize the 
#' distribution of C-scores associated with sites annotated in a RiboClass object.


.plot_boxplot_sites <- function(matrix, values_to_plot, outlier, sort_by, horizontal) { # main function called in the Rmd
  site <- cscore <- NULL # variable init
  id_vars <- "site" # usfull for long format transformation
    
    # Site sorting according to 3 criteria (iqr, median or var) by default: median
    if (tolower(sort_by) == "iqr") { # convert the variable sort_by to a minuscule and verifies if sort_by = "iqr"
      matrix_melted <- get_IQR(matrix) # if TRUE, the get_IQR function is called
      
    } else if(tolower(sort_by) == "median") { # otherwise converts the variable sort-by to lowercase and checks if sort_by = "median"
      matrix_melted <- reshape2::melt(matrix, id.vars = id_vars, # Transformation of the table matrix into a long format id.vars = site
                                      value.name = values_to_plot) # Column that will store the values ("cscore")
      print(head(matrix_melted))
      
    } else if (tolower(sort_by) == "var")  { # otherwise if sort_by = var we sort by var 
      matrix_melted <- get_IQR(matrix, "var") # function call : get_IQR 
      
    } else {
      stop("Choose either \"median\", \"iqr\" or \"var\" for sort_by param.\n  \"", # Error message if none of the 3 is recognized 
           sort_by,"\" is not a valid argument.")
    }
  
  shape_outlier <- NA 
  if (outlier) # if outlier = TRUE
    shape_outlier <- 19 # represented by full circle (19)
  
  if ((tolower(sort_by) %in% c("iqr","var"))) { # si sort_by = iqr ou var
    method <- ifelse(tolower(sort_by) =="iqr","IQR","variance") # dynamic title
   
  p <- ggplot(matrix_melted, aes(x = site, y = cscore)) + # plot creation
    geom_boxplot() +
    theme_bw() +
    ggtitle("Boxplot of C-score values across the samples") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust = 1),
          legend.position = "top",
          axis.title=element_text(size=14),
          axis.text.y = element_text(size=12)) +
    labs(x = paste0("rRNA 2'Ome sites (sorted by decreasing ",method,")"),
         y = "C-score",
         subtitle = paste0(
           length(unique(matrix_melted$sample)), " samples and ", 
           length(unique(matrix_melted$site)), " positions"
         )
    )
  return(p)
  
  } else { # otherwise (if sort_by = median)
    matrix_melted <- reshape2::melt(matrix, id.vars = id_vars, # sort site by median with stats::reorder()
                                    value.name = values_to_plot)
    
  p <- ggplot2::ggplot(matrix_melted,
                       ggplot2::aes(x = stats::reorder(site, # Sorts the site values based on a metric (metric)
                                                       !!rlang::sym(values_to_plot), # dynamically retrieves the corresponding column to values_to_plot ("cscore", "median", "iqr", etc.)
                                                       na.rm = TRUE),
                                    y = !!rlang::sym(values_to_plot))) +
    ggplot2::geom_boxplot(outlier.shape = shape_outlier) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5,
                                                       hjust = 1)) +
    ggplot2::xlab("rRNA 2'Ome sites (sorted by median)") +
    ggplot2::ylab("C-score")
  }
  
  if (horizontal) # if horizontal = TRUE
    p <- p + ggplot2::coord_flip() # invert axis
  
  return(p)
}

#' (internal) plot boxplot for a given matrix of values
#'
#' @param matrix Sites x Samples C-score/count matrix (output of extract_data()).
#' @param metadata Metadata of samples in matrix.
#' @param color_col Name of the column in the metadata used for coloring samples.
#' @param outlier Show boxplot outlier values.
#' @param values_col_name Name of the column containing the value (either count
#' or cscore).
#' @param horizontal Show boxplot horizontally.
#'
#' @return ggplot boxplot
#' @keywords internal

.plot_boxplot_samples <- function(matrix, values_col_name, metadata,
                                  color_col = NA, outlier= TRUE, horizontal=TRUE) {
  Sample <- NULL
  id_vars <- "Sample"
  matrix <- log10(matrix) # conversion of c-score in log10
  matrix_inv <- as.data.frame(t(matrix)) # transposition of the matrix
  # matrix_inv <- tibble::rownames_to_column(matrix_inv, "Sample")
  matrix_inv["Sample"] <- rownames(matrix_inv)
  
  # verify if color_col is given
  if (!is.na(color_col)) {
    matrix_inv <- cbind(matrix_inv, metadata[color_col])
    id_vars <- c(color_col, "Sample")
  }else{
    # If color_col is missing, create a "qc" column to identify the outliers
    matrix_inv <- cbind(matrix_inv, qc=0)
    matrix_inv$qc[apply(matrix,2,median,"na.rm"=TRUE)<=2] <- 1 
    id_vars <- c("qc", "Sample")
    
  }
  # Transformation in format "melted" for ggplot
  matrix_melted <- reshape2::melt(matrix_inv, id.vars = id_vars, #creation of a matrix with "qc" columns (=0 or 1 if outlier) "Sample", "variable" et "count"
                                  value.name = values_col_name)
  # Transform "Sample" in factor to control order
  matrix_melted[["Sample"]] <- factor(matrix_melted[["Sample"]],
                                      levels = unique(matrix_melted[["Sample"]]))
  shape_outlier <- NA

  if (outlier)
    shape_outlier <- 19
  
  # Initialize the plot with or without the defined color
  if (is.na(color_col)) {   # if no color specified
    p <- ggplot2::ggplot(matrix_melted,
                         ggplot2::aes(x = Sample, y = !!sym(values_col_name), # create the qc column to map colors
                                      fill = factor(qc))) +
      ggplot2::geom_boxplot(outlier.shape = shape_outlier) +
      ggplot2::scale_fill_manual(values = c("0" = "white", "1" = "red")) +
      ggplot2::theme(legend.position = "none")
  } else {
    p <- ggplot2::ggplot(matrix_melted,
                         ggplot2::aes(x = Sample, y = !!sym(values_col_name), fill = !!sym(color_col))) +
      ggplot2::geom_boxplot(outlier.shape = shape_outlier)
  }
  #add boxplot
  p <- p + ggplot2::geom_boxplot(outlier.shape = shape_outlier) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5,
                                                       hjust = 1)) +
    ylab("log10(total end read count)")

  # add horizontal line
  p <- p + ggplot2::geom_hline(yintercept = 2, colour = "blue")

  # horizontal orientation (optional)
  if (horizontal)
    p <- p + ggplot2::coord_flip()

  return(p)
}