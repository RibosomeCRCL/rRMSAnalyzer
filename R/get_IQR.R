#' Compute IQR or variance for 2'Ome sites
#'
#' @param df.meth.scores a dataframe 
#' @param order athe type of variance to be computed. Should be var or iqr.
#' @return a ggplot data frame with Cscore, Site, IQR or var value and sample ID
#' @export
#'
#' @examples get_IQR(aggregate_samples_by_col(ribo[["counts"]],"cscore_median",position_to_rownames = T))


get_IQR <- function(df.meth.scores = NULL, order = "IQR"){
  
  
  #the df should have in rows 2'OMe sites and in columns the samples
  df.meth.scores <- as.data.frame(df.meth.scores)
  

  if (tolower(order) == "iqr") {
    df.meth.scores$iqr <- apply(df.meth.scores, 1, IQR) # compute the IQR
    df.meth.scores$sites.id <- rownames(df.meth.scores)
    
    df.meth.scores <- tidyr::gather(df.meth.scores, "sampleID", "Cscore", -sites.id, -iqr) # Gather for site id and IQR
    
    df.meth.scores <- df.meth.scores[order(df.meth.scores$iqr),] # order by IQR
    
    df.meth.scores$sites.id <- factor(df.meth.scores$sites.id, levels = unique(df.meth.scores$sites.id)) # refactor the site.id
  }
  
  if (tolower(order) == "var") {
    
    df.meth.scores$var <- apply(df.meth.scores, 1, var) # compute the var
    df.meth.scores$sites.id <- rownames(df.meth.scores)
    
    df.meth.scores <- tidyr::gather(df.meth.scores, "sampleID", "Cscore", -sites.id, -var) # Gather for site id and var
    
    df.meth.scores <- df.meth.scores[order(df.meth.scores$var),] # order by var
    
    df.meth.scores$sites.id <- factor(df.meth.scores$sites.id, levels = unique(df.meth.scores$sites.id)) # refactor the site.id

  }
  
  return(df.meth.scores)
  
}


#' Plot IQR or variance per site ID / Plot boxplot per site id 
#'
#' @param ribo a ribo class
#' @param plot boxplot or IQR (by defalut)
#' @param variance IQR or var for variance
#' 
#' @return a ggplot
#' @export
#'
#' @examples plot_sites_by_IQR(ribo = ribo.meth, plot = "boxplot")
#' plot_sites_by_IQR(ribo = ribo.meth)


plot_sites_by_IQR <- function(ribo = NULL, plot = "IQR", variance = "IQR") {
  
  Cscore_matrix <- aggregate_samples_by_col(ribo[["counts"]],"cscore_median",position_to_rownames = T)
  IQR_order_df <- get_IQR(Cscore_matrix, order = variance)
  
  
  if (tolower(plot) == "iqr") {
  
  p1 <- ggplot(IQR_order_df, aes(x = sites.id, y = iqr)) +
    geom_point() +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, size = 9, hjust = 1), 
          legend.position = "top", 
          axis.title=element_text(size=14), 
          axis.text.y = element_text(size=12)) + 
    labs(x = "rRNA 2'Ome sites",
         y = "Interquartile Range (IQR)", size = 20)
    
  }
  
  if (tolower(plot) == "boxplot") {
    
    p1 <- ggplot(IQR_order_df, aes(x = sites.id, y = Cscore)) +
      geom_boxplot() +
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 45, size = 9, hjust = 1), 
            legend.position = "top", 
            axis.title=element_text(size=14), 
            axis.text.y = element_text(size=12)) + 
      labs(x = "rRNA 2'Ome sites",
           y = "Cscore", size = 20)
    
  }
  
  return(p1)
}
