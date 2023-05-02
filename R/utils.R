check_metadata <- function(ribo,metadata_name) {
  # Get metadata columns that do not exist in ribo's metadata.
  unmatched_elts <- metadata_name[
    which(!(metadata_name %in% colnames(ribo[["metadata"]])))
    ]
  
  if(length(unmatched_elts) > 0) {
    stop("The following metadata columns you have specified do not exist : \n- ",
         paste(unmatched_elts,collapse = "\n- "),
         "\n\n Available columns : ",
       paste(colnames(ribo[["metadata"]]),collapse = ", "))
    
  }
}

# Adapted from: https://martinctc.github.io/blog/vignette-downloadable-tables-in-rmarkdown-with-the-dt-package/
create_dt <- function(x){
  DT::datatable(x,
                rownames = FALSE,
                extensions = 'Buttons',
                options = list(dom = 'Blfrtip',
                               buttons = list('copy', 'csv', list(
                                 extend='excel',
                                 text="Excel"), 'pdf')
                               ))
}
