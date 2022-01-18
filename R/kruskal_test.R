# Perform a Kruskal-Wallis rank sum test on C-score
# Takes in input a C-score matrix 
# The metadata
# Test to be done in rows or in columns

kruskal.test.on.cscores <- function(cscore.matrix = NULL, metadata = NULL, column.or.row = "row", order.by.col = NULL, adjust.pvalues.method = "fdr", factor.column = NULL) {
  # Quality check for each parameter
  cscore.matrix <- as.data.frame(cscore.matrix)
  if(is.null(cscore.matrix)) {stop("Missing cscore matrix. Please specify a cscore.matrix")}
  if(is.null(metadata)) {stop("Missing metadata. Please specify a metadata data frame")}
  if(is.null(order.by.col)) {stop("Missing parameter. Please specify a column to order cscore matrix")}
  
  if (tolower(column.or.row) == "row") {
    cscore.matrix <- cscore.matrix[complete.cases(cscore.matrix), match(metadata[,order.by.col], colnames(cscore.matrix))] # order column cscores as metadata
    kruskal.test.pvalues <- apply(cscore.matrix, 1, function(x) {kruskal.test(x ~ metadata[,factor.column])$p.value}) # apply to each row the test and extract p.value
  }
  
  if (tolower(column.or.row) == "column") {
    cscore.matrix <- cscore.matrix[match(metadata[,order.by.col], rownames(cscore.matrix)),] # order rows as metadata
    kruskal.test.pvalues <- apply(cscore.matrix, 2, function(x) {kruskal.test(x ~ metadata[,factor.column])$p.value}) # apply to each column the test and extract p.val
  }
  
  # a resume data frame
  # siteID
  # pvalue of the test
  # adjusted p.value with the same method as the parameter. For more adjsutment methods, please see ?p.adjust in R console
  df.kruskal.pvalues <- data.frame(siteID = names(kruskal.test.pvalues), p.val = kruskal.test.pvalues, 
                                   p.adj = p.adjust(kruskal.test.pvalues, method = adjust.pvalues.method))
  
  # retur the data frame with p values
  return(df.kruskal.pvalues)
  
}



mean_max_difference <- function(cscore.matrix = NULL, metadata = NULL, order.by.col = NULL, group.by.col = NULL) {
  # cscore matrix => 2'Ome sites in column and samples in row
  if(is.null(cscore.matrix)) {stop("Missing cscore matrix. Please specify a cscore.matrix")}
  if(is.null(metadata)) {stop("Missing metadata. Please specify a metadata data frame")}
  if(is.null(order.by.col)) {stop("Missing parameter. Please specify a column to order cscore matrix")}
  if(is.null(group.by.col)) {stop("Missing parameter. Please specify a column to group cscores")}
  
  cscore.matrix <- as.data.frame(cscore.matrix[match(metadata[,order.by.col], rownames(cscore.matrix)),]) #order samples according to order.by.col parameter
  
  df.mean.each.group <- aggregate(cscore.matrix, list(metadata[,group.by.col]), mean) # aggregate each site, compute the mean according to the group.by.col 
  
  df.min.max <- apply(df.mean.each.group[,-1], 2, function(x) {diff(range(x))}) # compute the difference between the max - min group 
  
  df.min.max <- data.frame(siteID = names(df.min.max), mean_max_min_difference = df.min.max) # df summary
  return(df.min.max)
}


wrapper.kruskal.test <- function(cscore.matrix = NULL, metadata = NULL, column.or.row = "row", order.by.col = NULL, adjust.pvalues.method = "fdr", factor.column = NULL) {
  
  df.pval <- kruskal.test.on.cscores(cscore.matrix = cscore.matrix, metadata = metadata, order.by.col = order.by.col , factor.column = factor.column, column.or.row = column.or.row)
  
  df.min.max <- mean_max_difference(cscore.matrix = cscore.matrix, metadata = metadata, order.by.col = order.by.col, group.by.col = factor.column)
  
  df.final <- merge(df.pval, df.min.max, by = "siteID")
  return(df.final)
  
}
