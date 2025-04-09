#' Return the most/less variant sites of a dataframe or a RiboClass
#'
#' @param ribo a RiboClass object.
#' @param n Number of top sites to return.
#' @param type_of_variant Either "less" or "most", to select the top n less
#' variant sites or the top n most variant sites respectively.
#' @param only_annotated (RiboClass only) Check variability only among
#' annotated sites. Ignored when df is a dataframe.
#'
#' @return A dataframe with the n most/less variant sites
#' @export
#'
#' @examples
#' data("ribo_toy")
#' data("human_methylated")
#' ribo_toy <- rename_rna(ribo_toy)
#' ribo_toy <- annotate_site(ribo_toy,human_methylated)
#' get_variant_sites(ribo = ribo_toy)
get_variant_sites <- function(ribo, n = NULL, type_of_variant = "most",
                              only_annotated = TRUE) {
    site <- NULL
    
    # Use n defined in plot_site_by_iqr.R
    if (exists("n_sites_variant", envir = .GlobalEnv)) {
      n <- get("n_sites_variant", envir = .GlobalEnv)
    } else {
      n <- 20  #Default values
    }
    
    df <- extract_data(ribo,only_annotated = only_annotated,
                       position_to_rownames = TRUE)
    
    df <- as.data.frame(df)
    
    if (tolower(type_of_variant) == "most") {
      most_variant <- TRUE
    } else if (tolower(type_of_variant) == "less") {
      most_variant <- FALSE
    } else {
      stop("Unrecognized type of variants. Options are \"most\" or \"less\"")
    }
    
      var <- apply(df, 1, var)
      df <- df[order(var, decreasing = most_variant)[1:n],]
      df["site"] <- rownames(df)
      df <- dplyr::relocate(df,site)
    return(df)
    
}