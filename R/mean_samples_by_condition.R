#' Regroup samples by condition and calculate mean for each condition
#' 
#' @description An helper function that will give the mean of the cscore or count of all samples by condition, for each position. The standard deviation is also given.
#' this can be used to create boxplot with ggplot.
#'
#' @param ribo a RiboClass object
#' @param metadata_condition name or position of the column __in metadata__ containing the condition
#' @param value name or position of the column containing the values on which mean by condition is calculated.
#' @param only_annotated use annotation site name instead of default
#' 
#' @importFrom dplyr %>%
#' @importFrom rlang sym
#' @return a dataframe with the mean for each condition for a selected value
#' @export
#' @md
#' @examples
#' data("ribo_toy")
#' mean_df <- mean_samples_by_conditon(ribo_toy,"count","condition")
mean_samples_by_conditon <- function(ribo,value, metadata_condition, only_annotated = F) {
  
  named_position <- NULL # NSE fix
  ribo_list <- ribo[["data"]]
  ribo_names <- names(ribo_list)
  ribo_list_named <- lapply(ribo_names, function(x){
    ribo_list[[x]]["sample"] <- x
    return(ribo_list[[x]])
  })
  
  ribo_concat <- dplyr::bind_rows(ribo_list_named)
  
  metadata <- ribo[["metadata"]]
  ribo_concat[metadata_condition] <- metadata[,metadata_condition][match(ribo_concat[,"sample"], metadata[,"samplename"])]
  ribo_condition <- ribo_concat %>% dplyr::group_by(named_position, !!sym(metadata_condition)) %>% dplyr::summarise(mean = mean(!!sym(value)), sd = stats::sd(!!sym(value)))
  return(ribo_condition)
}