check_metadata <- function(ribo,metadata_name) {
  
  if(!(metadata_name %in% colnames(ribo[["metadata"]]))) {
    stop("The metadata column specified (",metadata_name,") does not exist.\n Available columns : ",
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
