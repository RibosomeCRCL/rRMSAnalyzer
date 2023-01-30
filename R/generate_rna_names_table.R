# (internal) Generate a table with former and current rna names for a given RiboClass
.generate_rna_names_table <- function(count_df) {
  rna_counts <- as.data.frame(sort(table(count_df[, 1])))
  rna_names_df <- data.frame(original_name = rna_counts[[1]], current_name =rna_counts[[1]])
  return(rna_names_df)
}