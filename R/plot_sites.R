#' Return a list of differential sites
#' 
#' @inherit plot_diff_sites details
#' @param df_of_kruskal Output of kruskal_test_on_cscores().
#' @param p_cutoff Cutoff below which the kruskal-wallis test is considered
#' significant for a given site.
#' @param cscore_cutoff Cutoff above which the max-min c-score range between
#' conditions' mean is considered significant.
#'
#' @return a vector of sites that are differentials
#' @keywords internal

select_most_differential_sites <- function(df_of_kruskal = NULL, p_cutoff = 1e-02, cscore_cutoff = 0.05){
  most_differential <- which(df_of_kruskal$p.adj < p_cutoff &
                                     df_of_kruskal$mean_max_min_difference > cscore_cutoff )
  most_differential_sites <- df_of_kruskal$site[most_differential]
  message(length(most_differential_sites), " significant sites found !")
  return(most_differential_sites)
}

#' (Internal) Display differential sites
#' 
#' Using dataframes obtained from
#' 
#' @param df_of_Cscores Dataframe of positions x samples for c-score.
#'  Output of extract_data().
#' @param df_of_kruskal Output of kruskal_test_on_cscores()
#' @param most_differential_sites A list of most differential sites.
#'  Output of select_most_differential_sites function.
#' @import ggplot2
#' @return a ggplot object
#' @keywords internal
#'
plot_most_differential_sites <- function(df_of_Cscores = NULL,
                                         df_of_kruskal = NULL,
                                         most_differential_sites = NULL) {
   group.id <- Cscore <- p.adj <- NULL

   df_of_Cscores <- tidyr::pivot_longer(df_of_Cscores,
                                        cols = c(2:ncol(df_of_Cscores)))
   
   colnames(df_of_Cscores) <- c("site","group.id","Cscore")
   df_of_Cscores$group.id <- gsub("\\..*","",df_of_Cscores$group.id)
   
   # df_of_cscores has the following 
   # group.id          siteID    Cscore
   # 1 Astrocytoma (A) 18S_Am99 0.9685661
   # 2 Astrocytoma (A) 18S_Am99 0.9599728
   # 3 Astrocytoma (A) 18S_Am99 0.9664454
   # 4 Astrocytoma (A) 18S_Am99 0.9743840
   # 5 Astrocytoma (A) 18S_Am99 0.9460999
   # 6 Astrocytoma (A) 18S_Am99 0.9656156
   
   cscore_diff_sites <- df_of_Cscores[which(df_of_Cscores$site %in% most_differential_sites),]
     p1 <- ggplot2::ggplot(cscore_diff_sites,ggplot2::aes(x = group.id, y = Cscore, fill = group.id)) + 
    geom_boxplot() +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x=element_blank(),
          legend.position = "top",
          legend.text = element_text(size = 12),
          axis.title.y = element_text(size = 16),
          text = element_text(size = 14)) +
    ylim(0,1) +
    facet_wrap(~site, nrow = 2) +
       geom_text(
         data = df_of_kruskal[which(df_of_kruskal$site %in% most_differential_sites), ],
         aes(
           x = -Inf,
           y = -Inf,
           label = paste(
             "p = ",
             scales::scientific(p.val, digits = 3) ,
             " \np.adj = ",
             scales::scientific(p.adj, digits = 3),
             sep = ""
           )
         ),
         hjust = -0.1,
         vjust = -1,
         inherit.aes = FALSE,
         size = 4
       ) +
       labs(x = "",
            y = "C-score",
            fill = "")

  return(p1)
}

#' Plot sites showing statistical difference between conditions
#' 
#' Display a boxplot between condition for each site where there is a sufficent difference between conditions.
#' The Kruskal-Wallis p-value is displayed at the bottom left of each boxplot.
#'
#' @details 
#' To be considered as differential, a site must follow two conditions :
#'  
#'   - Have a significant p-value on a Kruskal-wallis
#'    test on c-score between conditions.
#'   - Have a C-score range between conditions (max median - min mea) 
#'   above a certain cutoff.
#'
#' Both the p-value cutoff and c-score range cutoff can be changed with p_cutoff and cscore_cutoff parameters respectively.
#' @md
#' @param ribo A RiboClass object.
#' @param factor_column Metadata column used to group samples by.
#' @param p_cutoff Cutoff for the adjusted p-value of the kruskal-wallis test.
#' @param cscore_cutoff Cutoff above which the max-min c-score range between conditions' mean is considered significant.
#' @param adjust_pvalues_method Method used to adjust p-value (one of p.adjust.methods)
#' @param object_only Return the results of the kruskal-wallis and C-score mean range in a dataframe directly.
#' @return a ggplot object.
#' @export
#'
#' @examples
#' data("ribo_toy")
#' data("human_methylated")
#' ribo_toy <- rename_rna(ribo_toy)
#' ribo_toy <- annotate_site(ribo_toy,human_methylated)
#' plot_diff_sites(ribo_toy,"condition", p_cutoff=0.1)
plot_diff_sites <- function(ribo, factor_column,
                            p_cutoff = 0.05,
                            cscore_cutoff = 0.05,
                            adjust_pvalues_method = "fdr",
                            object_only = FALSE) {
  
  site <- NULL
  
  check_metadata(ribo,factor_column)
  ribo_matrix <- extract_data(ribo, only_annotated = TRUE, position_to_rownames = TRUE)
  kruskal_df <- wrapper_kruskal_test(ribo, factor_column = factor_column)  
  
 if(object_only) return(kruskal_df)
  
  #we replace sample names with their group name (the values in factor_column) 
  new_colnames <- lapply(colnames(ribo_matrix), function(x) {
    x <- ribo[["metadata"]][[factor_column]][ribo$metadata$samplename == x]
    return(x)
  })

  colnames(ribo_matrix) <- new_colnames
  ribo_matrix_rn <- ribo_matrix
  # Because ggplot is not a fan of rownames, we transfer the latter to its own column
  ribo_matrix_rn["site"] <- rownames(ribo_matrix)
  # Put sites column at first position.
  ribo_matrix_rn <- dplyr::relocate(ribo_matrix_rn,site)
  
  most_signi <- select_most_differential_sites(kruskal_df,p_cutoff = p_cutoff, cscore_cutoff = cscore_cutoff)
  if (length(most_signi) == 0) {
    return(ggplot() + 
      annotate("text", x = 4, y = 25, size=8,
               label = "No differential site found !") + 
      theme_void())
  }
  return(plot_most_differential_sites(df_of_Cscores = ribo_matrix_rn,
                                     df_of_kruskal = kruskal_df,
                                     most_differential_sites = most_signi))
  
}

