# Perform a Kruskal-Wallis rank sum test on C-score
# Takes in input a C-score matrix 
# The metadata
# Test to be done in rows or in columns

kruskal_test_on_cscores <- function(cscore_matrix = NULL, metadata = NULL, order_by_col = NULL, adjust_pvalues_method = "fdr", factor_column = NULL) {
  # Quality check for each parameter
  cscore_matrix <- as.data.frame(cscore_matrix)
  if(is.null(cscore_matrix)) {stop("Missing cscore matrix. Please specify a cscore_matrix")}
  if(is.null(metadata)) {stop("Missing metadata. Please specify a metadata data frame")}
  if(is.null(order_by_col)) {stop("Missing parameter. Please specify a column to order cscore matrix")}
  
  

    cscore_matrix <- cscore_matrix[complete.cases(cscore_matrix), match(metadata[,order_by_col], colnames(cscore_matrix))] # order column cscores as metadata
    kruskal_test_pvalues <- apply(cscore_matrix, 1, function(x) {kruskal.test(x ~ metadata[,factor_column])$p.value}) # apply to each row the test and extract p.value
  
  

  
  # a resume data frame
  # siteID
  # pvalue of the test
  # adjusted p.value with the same method as the parameter. For more adjsutment methods, please see ?p.adjust in R console
  df_kruskal_pvalues <- data.frame(siteID = names(kruskal_test_pvalues), p.val = kruskal_test_pvalues, 
                                   p.adj = p.adjust(kruskal_test_pvalues, method = adjust_pvalues_method))
  
  # retur the data frame with p values
  return(df_kruskal_pvalues)
  
}



mean_max_difference <- function(cscore_matrix = NULL, metadata = NULL, order_by_col = NULL, group.by.col = NULL) {
  # cscore matrix => 2'Ome sites in column and samples in row
  if(is.null(cscore_matrix)) {stop("Missing cscore matrix. Please specify a cscore_matrix")}
  if(is.null(metadata)) {stop("Missing metadata. Please specify a metadata data frame")}
  if(is.null(order_by_col)) {stop("Missing parameter. Please specify a column to order cscore matrix")}
  if(is.null(group.by.col)) {stop("Missing parameter. Please specify a column to group cscores")}
  
  cscore_matrix <- t(cscore_matrix)
  cscore_matrix <- as.data.frame(cscore_matrix[match(metadata[,order_by_col], rownames(cscore_matrix)),]) #order samples according to order_by_col parameter
  
  df.mean.each.group <- aggregate(cscore_matrix, list(metadata[,group.by.col]), mean) # aggregate each site, compute the mean according to the group.by.col 
  
  df.min.max <- apply(df.mean.each.group[,-1], 2, function(x) {diff(range(x))}) # compute the difference between the max - min group 
  
  df.min.max <- data.frame(siteID = names(df.min.max), mean_max_min_difference = df.min.max) # df summary
  return(df.min.max)
}

#' Title
#'
#' @param cscore_matrix 
#' @param metadata 
#' @param column_or_row 
#' @param order_by_col 
#' @param adjust_pvalues_method 
#' @param factor_column 
#'
#' @return
#' @export
#'
#' @examples
wrapper.kruskal.test <- function(ribo = NULL, order_by_col = NULL, adjust_pvalues_method = "fdr", factor_column = NULL) {
  
  metadata <- ribo[["metadata"]]
  cscore_matrix <- aggregate_samples_by_col(ribo$counts,"cscore_median",position_to_rownames = T) 
  df.pval <- kruskal_test_on_cscores(cscore_matrix = cscore_matrix, metadata = metadata, order_by_col = order_by_col , factor_column = factor_column)
  
  df.min.max <- mean_max_difference(cscore_matrix = cscore_matrix, metadata = metadata, order_by_col = order_by_col, group.by.col = factor_column)
  
  df.final <- merge(df.pval, df.min.max, by = "siteID")
  return(df.final)
  
}
