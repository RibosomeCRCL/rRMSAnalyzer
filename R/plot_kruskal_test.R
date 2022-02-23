#' Title
#'
#' @param df_of_kruskal 
#' @param pCutoff 
#' @param CscoreCutoff 
#'
#' @return
#' @export
#'
#' @examples
select_most_significant_sites <- function(df_of_kruskal = NULL, pCutoff = 1e-02, CscoreCutoff = 0.05){
  most_significant_sites <- df_of_kruskal$siteID[which(df_of_kruskal$p.adj < pCutoff & df_of_kruskal$mean_max_min_difference > CscoreCutoff )]
  return(most_significant_sites)
}

#' Title
#'
#' @param df_of_Cscores 
#' @param df_of_kruskal 
#' @param most_significant_sites 
#'
#' @return
#' @export
#'
#' @examples
plot_most_significant_sites <- function(df_of_Cscores = NULL, df_of_kruskal = NULL, most_significant_sites = most_significant_sites) {


   df_of_Cscores <- tidyr::pivot_longer(df_of_Cscores,cols = c(2:ncol(df_of_Cscores))) 
   colnames(df_of_Cscores) <- c("siteID","group.id","Cscore")
   df_of_Cscores$group.id <- gsub("\\..*","",df_of_Cscores$group.id)
   
   # df_of_cscores has the following 
   # group.id          siteID    Cscore
   # 1 Astrocytoma (A) 18S_Am99 0.9685661
   # 2 Astrocytoma (A) 18S_Am99 0.9599728
   # 3 Astrocytoma (A) 18S_Am99 0.9664454
   # 4 Astrocytoma (A) 18S_Am99 0.9743840
   # 5 Astrocytoma (A) 18S_Am99 0.9460999
   # 6 Astrocytoma (A) 18S_Am99 0.9656156

     p1 <- ggplot(df_of_Cscores[which(df_of_Cscores$siteID %in% most_significant_sites),] , aes(x = group.id, y = Cscore, fill = group.id)) + 
    geom_boxplot() +  
    theme_bw() + 
    theme(axis.text.x = element_blank(), 
          axis.ticks.x=element_blank(), 
          legend.position = "top",
          legend.text = element_text(size = 12), 
          axis.title.y = element_text(size = 16), 
          text = element_text(size = 14)) + 
    ylim(0,1) + 
    facet_wrap(~siteID, nrow = 2) +
    geom_text(data = df_of_kruskal[which(df_of_kruskal$siteID %in% most_significant_sites),], 
              aes(x = -Inf, y = -Inf, label = paste("p = ",scientific(p.adj, digits = 3), sep = "")), 
              hjust=-0.1, 
              vjust=-1, 
              inherit.aes = FALSE, 
              size=4) + 
    labs(x = "",
         y = "C-score",
         fill = "")
  
  return(p1)
}

#' Title
#'
#' @param ribo 
#' @param col_to_test 
#' @param values_column 
#' @param factor_column 
#'
#' @return
#' @export
#'
#' @examples
plot_test_wrapper <- function(ribo, values_column , factor_column, p_cutoff = 1e-02, cscore_cutoff = 0.05) {
  
  
  ribo_matrix <- aggregate_samples_by_col(ribo[["counts"]],values_column,position_to_rownames = T)
  kruskal_df = wrapper.kruskal.test(ribo,order_by_col = "samplename",factor_column = factor_column)

  #we replace sample names with their group name
  new_colnames <- lapply(colnames(ribo_matrix), function(x) {
    x <- ribo[["metadata"]][[factor_column]][ribo$metadata$samplename == x]
    return(x)
  })
  
  colnames(ribo_matrix) <- new_colnames
  ribo_matrix_rn <- ribo_matrix
  # Because GGplot is not a fan of rownames, we transfer the latter to its own column
  ribo_matrix_rn["named_position"] <- rownames(ribo_matrix)
  # Put named_position at first position.
  ribo_matrix_rn <- dplyr::select(ribo_matrix_rn,"named_position",tidyselect::everything())

  #TODO : Add metadata column directly to ribo_matrix

  
  
  most_signi <- select_most_significant_sites(kruskal_df,pCutoff = p_cutoff, CscoreCutoff = cscore_cutoff)
  
  plot_most_significant_sites(df_of_Cscores = ribo_matrix_rn,df_of_kruskal = kruskal_df,most_significant_sites = most_signi)
  
}
