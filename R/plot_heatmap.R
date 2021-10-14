#' plot heatmap for a riboclass object
#' @description This easy function will let you display an heatmap for any given column (count or c-score). You add an additionnal layer of information with metadata columns.
#' @param ribo a riboclass object
#' @param col_to_plot column containing the values for the heatmap
#' @param order_by_col column that contains your samplenames. Makes the link between your data and metadata
#' @param cols_for_annotation (optional) metadata columns to show on the heatmap
#' @param most_variant boolean to select the most variants
#'
#' @return
#' @export
#'
#' @examples plot_heatmap(ribo_with_score, "ScoreC.Median.net", cols_for_annotation = c("group", "group_num"), order_by_col = "samplename")
plot_heatmap <- function(ribo, col_to_plot, cols_for_annotation, order_by_col, most_variant = F) {
  matrix <- aggregate_samples_by_col(ribo[["raw_counts"]], col_to_plot, position_to_rownames = T)
  .plot_heatmap(matrix, ribo[["metadata"]], order_by_col, cols_for_annotation = cols_for_annotation, most_variant = most_variant)
}

#' Title
#'
#' @param cscore_matrix
#' @param annotation_samples
#' @param order_by_col
#' @param cols_for_annotation
#' @param most_variant
#'
#' @return
#' @export
#'
#' @examples
.plot_heatmap <- function(cscore_matrix = NULL, annotation_samples = NULL, cols_for_annotation = NULL, order_by_col = NULL, most_variant = F) {
  annotation_samples <- annotation_samples[, c(order_by_col, cols_for_annotation)]
  # rownames of annotation_samples
  rownames(annotation_samples) <- annotation_samples[, order_by_col]


  # all character columns to factor columns
  annotation_samples[sapply(annotation_samples, is.character)] <- lapply(
    annotation_samples[sapply(annotation_samples, is.character)],
    as.factor
  )

  annotation_samples_1 <- annotation_samples[cols_for_annotation]

  if (most_variant) {
    cscore_matrix <- get_most_or_less_variant(cscore_matrix)
  }


  htmap <- pheatmap::pheatmap(cscore_matrix[complete.cases(cscore_matrix), match(annotation_samples[, order_by_col], colnames(cscore_matrix))],
    clustering_method = "ward.D2",
    clustering_distance_cols = "manhattan",
    clustering_distance_rows = "manhattan",
    #    color = viridis::magma(100),
    cutree_cols = 2,
    cutree_rows = 4,
    annotation_col = data.frame(annotation_samples_1), main = "heatmap"
  )
  return(htmap)
}

# ==================================================

#' plot a correlation heatmap from a riboclass object
#' @description This heatmap is used to compare
#'
#' @param ribo a riboclass object
#' @param col_to_plot column containing the values to be used for the heatmap.
#' @param cols_for_annotation metadata columns for annotation
#' @param order_by_col column that contains your sample names. Makes the link between your data and metadata
#' @param use_triangle Show only a correlation triange without any more values
#'
#' @return
#' @export
#'
#' @examples
plot_heatmap_corr <- function(ribo, col_to_plot, cols_for_annotation, order_by_col, use_triangle = F) {
  matrix <- aggregate_samples_by_col(ribo[["raw_counts"]], col_to_plot, position_to_rownames = T)
  .plot_heatmap_corr(matrix, ribo[["metadata"]], order_by_col = order_by_col, cols_for_annotation = cols_for_annotation, use_triangle)
}
#' Internal function of plot_heatmap_corr.
#'
#' @param cscore_matrix
#' @param annotation_samples
#' @param cols_for_annotation
#' @param order_by_col
#' @param use_triangle
#'
#' @return
#' @export
#'
#' @examples
.plot_heatmap_corr <- function(cscore_matrix = NULL, annotation_samples = NULL, cols_for_annotation = NULL, order_by_col = NULL, use_triangle = F) {

  # rownames of annotation_samples
  rownames(annotation_samples) <- annotation_samples[, order_by_col]


  # all character columns to factor columns
  annotation_samples[sapply(annotation_samples, is.character)] <- lapply(
    annotation_samples[sapply(annotation_samples, is.character)],
    as.factor
  )
  corr_matrix <- 1 - cor(cscore_matrix)
  dist_cor <- as.dist(corr_matrix)
  annotation_samples_1 <- annotation_samples[cols_for_annotation]

  white_red <- colorRampPalette(c("white", "red"), interpolate = "linear")(100)

  if (use_triangle) {
    corr_matrix[lower.tri(corr_matrix)] <- NA

    htmap <- pheatmap::pheatmap(corr_matrix[, match(annotation_samples[, order_by_col], colnames(corr_matrix))],
      cluster_cols = F,
      cluster_rows = F,
      color = white_red,
      breaks = seq(0, 1, by = 0.01),
      annotation_col = data.frame(annotation_samples_1), main = "correlation heatmap"
    )
  } else {
    htmap <- pheatmap::pheatmap(corr_matrix[, match(annotation_samples[, order_by_col], colnames(corr_matrix))],
      clustering_method = "ward.D2",
      clustering_distance_cols = dist_cor,
      clustering_distance_rows = dist_cor,
      cutree_cols = 2,
      cutree_rows = 4,
      color = white_red,
      breaks = seq(0, 1, by = 0.01),
      annotation_col = data.frame(annotation_samples_1), main = "correlation heatmap"
    )
  }
  return(htmap)
}

get_most_or_less_variant <- function(df = NULL, column_or_row = "row", n = 20, type_of_variant = "most") {
  df <- as.data.frame(df)

  if (column_or_row == "row") {
    column_or_row <- 1
  }

  if (column_or_row == "column") {
    column_or_row <- 2
  }

  if (is.null(df)) {
    stop("No matix in entry")
  }

  var <- apply(df, column_or_row, var)
  if (tolower(type_of_variant) == "most") {
    get_most_variant <- TRUE
  } else if (tolower(type_of_variant) == "less") {
    get_most_variant <- FALSE
  } else {
    stop("Unrecognized type of variants. Options are \"most\" or \"less\"")
  }
  df[order(var, decreasing = get_most_variant)[1:n], ]
}
