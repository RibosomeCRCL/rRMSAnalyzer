#' Title
#'
#' @param Ribo 
#'
#' @return
#' @export
#'
#' @examples
plot_count_profile_test <- function(ribo) {
  
  ribo_matrix <- aggregate_samples_by_col(sample_list = ribo[[1]],col_to_keep ="Count",position_to_rownames = F)
  ribo_matrix_melted <- reshape2::melt(ribo_matrix)
  ribo_matrix_melted$rna  <- stringr::str_extract(ribo_matrix_melted$named_position,
                                                                 paste0("^([^", "_","])+"))
  #ribo_matrix_melted$named_position <- 1:length(ribo_matrix_melted$named_position)
  QC_plot <- ggplot(data=ribo_matrix_melted, aes(y=value, x=named_position,group=variable,color=variable)) +
    facet_wrap(~rna,scales = "free") +
    theme_bw() +
    geom_line() +
    xlab("Genomic position")+
    ylab("Raw 5'/3'-end read count") +
    theme(axis.text.x = element_text(angle=90,size = 1))
  
  return(QC_plot)
  
  
}

#' Title
#'
#' @param ribo 
#' @param metadata_col 
#'
#' @return
#' @export
#'
#' @examples
plot_count_profile_by_cond <- function(ribo, metadata_col) {
  ribo_concat <- mean_samples_by_conditon(ribo,"Count",metadata_col )
  ribo_concat["rna"]  <- stringr::str_extract(ribo_concat$named_position,paste0("^([^", "_","])+"))
                                                  
  lineplot <- ggplot(ribo_concat, aes(x=named_position, y=mean, colour=!!sym(metadata_col), group = !!sym(metadata_col))) + 
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1, position = position_dodge(0.1)) +
    geom_line(position = position_dodge(0.1)) +
    geom_point(position = position_dodge(0.1), size=3, shape=21, fill="white") + 
    facet_wrap(~rna,scales = "free") +
    theme_bw() +
    geom_line() +
    xlab("Genomic position")+
    ylab("Raw 5'/3'-end read count") +
    theme(axis.text.x = element_text(angle=90,size = 1))
  
  
return(lineplot)
}