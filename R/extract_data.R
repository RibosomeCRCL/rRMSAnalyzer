#' Aggregate results into a single matrix
#'
#' For a given column in data, this function will generate a dataframe with all samples.
#' Exports all positions (if only_annotated is false) or only annotated positions
#' (if only_annotated is true).
#'
#' @param ribo A RiboClass object.
#' @param col Column in data you want extract data from (cscore or count).
#' @param position_to_rownames If true, position will be included as a rowname. 
#' They will in a new column otherwise.
#' @param only_annotated If true, return a dataframe with only annotated sites.
#' Return all sites otherwise.
#' @return A position x samples dataframe. Position can be all or only
#' annotated sites. A site column can be added if position_to_rownames = FALSE.
#' @export
#'
#' @examples
#' data('ribo_toy')
#' count_df <- extract_data(ribo_toy,'count')
extract_data <- function(ribo, col = "cscore",
                         position_to_rownames = FALSE, only_annotated = FALSE) {
  named_position <- NULL  # NSE fix
  # The rows of this matrix correspond to the positions on the rRNA
  col <- tolower(col)
  
  if (!(col %in% colnames(ribo[["data"]][[1]]))) {
    cli::cli_abort(c("Name supplied to {.var col} is not a column in the data",
                   "i" = "Available columns: {.val {colnames(ribo[[\"data\"]][[1]])}}",
                   "x" = "{.val {col}} is not a valid column"))
  }
  
  sample_list <- ribo[["data"]]
  sample_list_nm <- names(sample_list)
  if (only_annotated) {
    # if only_annotated is TRUE, then we use the "site" column for our
    # positions, as this column contains the annotated sites only.
    df <- sample_list[[1]]
    df_sites <- df[which(!is.na(df[, "site"])), ]
    position_list <- df_sites[, "site"]
    matrix_all <- data.frame(site = position_list)
    positions <- "site"
  } else {
    # When only_annotated is FALSE, we default to all positions in
    # "named_position".
    position_list <- sample_list[[1]][, "named_position"]
    matrix_all <- data.frame(named_position = position_list)
    positions <- "named_position"
    
  }
  
  for (sample_nm in sample_list_nm) {
    sample_df <- sample_list[[sample_nm]]
    if (only_annotated) {
      sample_df <- sample_df[which(!is.na(sample_df[, "site"])), ]
    }
    sample_df <- sample_df[, c(positions, col)]
    
    matrix_all <- dplyr::full_join(matrix_all, sample_df,
                                     by = positions)
    matrix_all_len <- length(names(matrix_all))
    names(matrix_all)[matrix_all_len] <- sample_nm
    matrix_all <- matrix_all[match(position_list,matrix_all[, positions]),]
    
  }
  if (position_to_rownames) {
    if (only_annotated) {
      col_name <- "site"
    } else {
      col_name <- "named_position"
    }
    row.names(matrix_all) <- matrix_all[, col_name]
    matrix_all[, col_name] <- NULL
  }
  return(matrix_all)
}