#' Compute IQR or variance for 2'Ome sites
#'
#' @param df a dataframe 
#' @param order athe type of variance to be computed. Should be var or iqr.
#' @return a ggplot data frame with Cscore, Site, IQR or var value and sample ID
#' @export
#'
#' @examples 
#' data("ribo_toy")
#' ribo_matrix <- extract_data(ribo_toy,position_to_rownames = TRUE)
#' ribo_matrix <- na.omit(ribo_matrix)
#' IQR_df <- get_IQR(ribo_matrix)
get_IQR <- function(df = NULL, order = "IQR"){
  
  IQR <- site <- iqr <- var <- NULL
  
  has_site_col <- FALSE # do the dataframe uses a site column ?
  
  #the df should have in rows 2'OMe sites and in columns the samples
  df <- as.data.frame(df)
  
  if ("site" %in% colnames(df)) {
    has_site_col <- TRUE
    rownames(df) <- df[["site"]]
    df["site"] <- list(NULL)
  }
  
  if (tolower(order) == "iqr") {
    df$iqr <- apply(df, 1, IQR) # compute the IQR
    df$site <- rownames(df)
    
    df <- tidyr::gather(df, "sample", "cscore", -site, -iqr) # Gather for site id and IQR
    
    df <- df[order(df$iqr, decreasing = TRUE),] # order by IQR
    
    df$site <- factor(df$site, levels = unique(df$site)) # refactor the site.id
  }
  
  if (tolower(order) == "var") {
    
    df$var <- apply(df, 1, var) # compute the var
    df$site <- rownames(df)
    
    df <- tidyr::gather(df, "sample", "cscore", -site, -var) # Gather for site id and var
    
    df <- df[order(df$var,decreasing = TRUE ),] # order by var
    
    df$site <- factor(df$site, levels = unique(df$site)) # refactor the site.id
    
  }
  
  return(df)
  
}