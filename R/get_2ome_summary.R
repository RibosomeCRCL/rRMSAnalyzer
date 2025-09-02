#' Create a metric summary table for 2ome_analysis report
#'
#' @param ribo a RiboClass object
#'
#' @returns a data frame
#' @export
#'
#' @examples
#' data("ribo_toy")
#' data("human_methylated")
#' ribo_toy <- remove_ribo_samples(ribo_toy,c("RNA1", "RNA2"))
#' ribo_toy <- rename_rna(ribo_toy)  
#' ribo_toy <- annotate_site(ribo_toy,
#'                                 annot = human_methylated,
#'                                 anno_rna = "rRNA",
#'                                 anno_pos = "Position",
#'                                 anno_value = "Nomenclature")
#' get_2ome_summary(ribo = ribo_toy)
get_2ome_summary <-function(ribo = ribo){
  #get IQR
  Cscore_matrix <- extract_data(ribo, col = "cscore", position_to_rownames = TRUE, only_annotated = TRUE)
  
  IQR_df <- rRMSAnalyzer::get_IQR(Cscore_matrix, order = "IQR")
  IQR_df <- IQR_df[order(IQR_df[,1], decreasing = TRUE), ]
  IQR_df$site <- factor(x = IQR_df$site, levels = unique(IQR_df$site))
  
  # IQR stats for each sites
  stats_IQR_df <- IQR_df %>%
    dplyr::group_by(site) %>%
    dplyr::summarise(
      IQR_median_all_samples = median(iqr))   # IQR median calculation
  
  # Extract data
  df_ribo <- ribo[[1]][[1]]
  variant_sites <- rRMSAnalyzer::get_variant_sites(ribo)[["site"]]
  df_ribo_annot2 <- df_ribo[which(!is.na(df_ribo["site"])),c("rna","rnapos","site","cscore")]
  
  df_ribo_annot2 <- df_ribo_annot2 %>%
    dplyr::rename(median_c_score_all_samples = cscore) %>% #replace c-score by median_c_score_all_samples
    dplyr::mutate(status_all_samples = ifelse(site %in% variant_sites, "variable", "stable")) %>% # add column status_all_samples to flag sites in red in "plot_site_by_IQR" 
    dplyr::select(-rna) #erase rna column 
  
  
  df_ribo_annot2 <- df_ribo_annot2 %>%
    dplyr::left_join(stats_IQR_df, by = "site")
  
  df_ribo_annot2 <- df_ribo_annot2 %>%
    dplyr::select(rnapos, site, status_all_samples, median_c_score_all_samples, IQR_median_all_samples)
}