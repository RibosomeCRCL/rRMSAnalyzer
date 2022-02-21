#' Transfrom ribo class data to ggplot data frame
#'
#' @param ribo 
#' @param metadata_columns The metadata columns to keep when transforming
#'
#' @return
#' @export
#'
#' @examples
#' transform.ribo(ribo = ribo, metadata_columns = c(5,6))
#' transform.ribo(ribo = ribo.meth, metadata_columns = c("SAMPLE..ORIGIN", "SPECIES"))




transform.ribo <- function(ribo, metadata_columns = NULL) {
  
  # First let's extract the data from ribo
  # TODO this function should be simplified, 
  # specially the values_column parameter which should be only "Cscore" or "Counts"
  values_column <- "cscore_median"
  df.matrix <- aggregate_samples_by_col(ribo[["counts"]],values_column,position_to_rownames = T)
  # Add a new columns that correspond to rownames
  # TODO make sure that aggregate_samples_by_col returns the good rownames (with the site Name)
  df.matrix$siteID <- rownames(df.matrix)
  # Transform the data with tidyr
  df.tranform <- tidyr::gather(df.matrix, "sampleID", "Cscore", -siteID)
  
  if (is.null(metadata_columns)) {
    # if no metadata, return the transformed data frame
    return(df.tranform)
  }
  else {
    # Check if the metadata_columns are numeric. 
    # If yes than get the column names of the selected columns
    if (is.numeric(metadata_columns)){
      metadata_columns <- colnames(ribo[["metadata"]][metadata_columns])
    }
    # Get the metadata
    # TODO samplename column should always exists in metadata
    df.meta <- ribo.meth[["metadata"]][c("samplename",metadata_columns)]
    # if metadata is empty, return an error
    if (nrow(df.meta) == 0) { stop("No metadata found") }
    # Merge the metadata with the transformed data frame
    df.tranform <- merge(df.tranform, df.meta, by.x = "sampleID", by.y = "samplename")
    return(df.tranform)
  }
}
