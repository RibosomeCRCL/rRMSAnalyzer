#' Create a outlier summary table (used in plot_heatmap_annotated)
#'
#' @param ribo a RiboClass object
#'
#' @returns a data frame
#' @export
#'
#' @examples
#' data("ribo_toy")
#' get_outliers(ribo = ribo_toy)
get_outliers <- function(ribo = ribo){

qcdata <- ribo[[2]]

## Create a df to save all the data
ribo_matrix <- rRMSAnalyzer::extract_data(ribo, "count", position_to_rownames = TRUE)
#calculate median for all the column in the ribo_matrix df
medcov <- apply(ribo_matrix,2,function(x) median(x,na.rm=TRUE)) 
# create a df like ribo_matrix and add median_coverage column
qcdata <- cbind(qcdata,median_coverage=medcov[rownames(qcdata)]) 
# add a quality column initialized
qcdata$coverage_quality <- "pass" 
# add a threshold for outliers
qcdata$coverage_quality[qcdata$median_coverage < 100] <- "warning" 

#plot object extraction
plot_rle <- rRMSAnalyzer::plot_rle(ribo,show_outlier = TRUE)
# extraction of the median absolute deviation
mad <- plot_rle$plot_env$mad  
# extraction of the median for all samples
rle_median <- as.vector(plot_rle$plot_env$rle_grouped$median)
# extraction of sample names
keys <- as.vector(plot_rle$plot_env$rle_grouped$key)
# creation of a data frame containing sample names and medians
outlier_table <- data.frame(
  key = keys,
  rle_median = rle_median)
# creation of the data frame outlier 
outlier <- data.frame()
# calculation of outliers 
for (i in 1:nrow(outlier_table)) {
  if(outlier_table$rle_median[i] < mad) {
    outlier <- rbind(outlier, outlier_table[i,])
  }
}

#incrementing the df qcdata with rle median
qcdata$rle_median <- outlier_table$rle_median[match(qcdata$samplename, outlier_table$key)]
# add a rle_median_quality column initialized at ok
qcdata$rle_median_quality <- "pass" 
# add a threshold for outliers
qcdata$rle_median_quality[qcdata$rle_median < mad] <- "warning" 

#incrementing qcdata
#calculate the number of fold the sample is an outlier
qcdata$total_outliers <- rowSums(qcdata[, c("coverage_quality", "rle_median_quality")] == "warning")

qcdata$total_outliers <- as.character(qcdata$total_outliers)

return(qcdata)
}