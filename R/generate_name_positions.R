#' Generate "named_position" column
#'
#' @param df dataframe with both rna and pos columns
#' @param rna name or position of rna column
#' @param pos name or position of position in rna column
#' @keywords internal
#' @return a dataframe
generate_name_positions <- function(df, rna, pos) {
  df["named_position"] <- paste(df[,rna], formatC(df[,pos],width = 4,flag = "0"), sep = "_")
  return(df)
}