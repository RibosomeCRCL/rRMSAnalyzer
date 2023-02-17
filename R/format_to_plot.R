#' Transform RiboClass data to a "ggplot-friendly" data frame.
#' 
#' Turn a riboClass into a dataframe that should be ready to be used in ggplot.
#' Metadata columns can be added if extra informations are needed.
#'
#' @param ribo A RiboClass object.
#' @param metadata_col Metadata columns to add. Must be in RiboClass' metadata.
#' @param only_annotated Keep only sites that have been annotated.
#' @return A ggplot-friendly dataframe with the following columns : 1) site; 2) sample; 3) cscore; ...) any metadata added
#' @export
#' 
format_to_plot <- function(ribo, metadata_col = NULL, only_annotated = FALSE) {
  site <- NULL #NSE fix
  # First let's extract the data from ribo
  values_column <- "cscore"
  df.matrix <- extract_data(ribo,values_column,position_to_rownames = T,
                            only_annotated = only_annotated)
  # Add a new columns that correspond to rownames
  df.matrix$site <- rownames(df.matrix)
  # Transform the data with tidyr
  df.tranform <- tidyr::gather(df.matrix, "sample", "cscore", -site)
  
  if (is.null(metadata_col)) {
    # if no metadata, return the transformed data frame
    return(df.tranform)
  }
  else {
    # Check if the metadata_columns are numeric. 
    # If true then get the column names of the selected columns
    if (is.numeric(metadata_col)){
      metadata_columns <- colnames(ribo[["metadata"]][metadata_col])
    }
    # Get the metadata
    df.meta <- ribo[["metadata"]][c("samplename",metadata_col)]
    # if metadata is empty, return an error
    if (nrow(df.meta) == 0) { stop("No metadata found") }
    # Merge the metadata with the transformed data frame
    df.tranform <- merge(df.tranform, df.meta, by.x = "sample", by.y = "samplename")
    return(df.tranform)
  }
}
