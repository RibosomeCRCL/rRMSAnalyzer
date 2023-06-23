# Generate the col argument for ComplexHeatmap::HeatmapAnnotation
generate_palette <- function(metadata,cols_to_use) {
  annot_list <- list()
  palettes_template <- list(
    c('#8dd3c7','#00bfff','#bebada','#f14292','#fdb462','#80b1d3','#b1de69','#fccde5','#d9d9d9'),
    c('#b2df8a','#ea1a8c','#a6cee3','#ffEf00','#1f78b4','#33a02c','#fddf6f','#cab2d6','#191970'),
    c('#fb9a99','#FFAA00','#AAFF00','#00FF00','#00FFAA','#00AAFF','#0000FF','#AA00FF','#FF00AA'))
  
  palettes <- palettes_template
  for(column in cols_to_use) {
    cond_names <- unique(metadata[[column]])
    if(is.numeric(cond_names)) col_is_numeric <- TRUE
    cond_names[which(is.na(cond_names))] <- "NA"
    annot <- c()
    
    if(any(col_is_numeric,length(cond_names) > 9)) {
      palette <- grDevices::hcl.colors(length(cond_names),"Light Grays")

    } else {
      palette <- palettes[[1]]
      
    }
    
    for(cond in cond_names) {
      annot[cond] <- palette[1]
      palette <- palette[-1]
      
    }
    annot_list[[column]] <- annot
    palettes <- palettes[-1]
    if(length(palettes) == 0) palettes <- palettes_template
  }
  return(annot_list)
}

