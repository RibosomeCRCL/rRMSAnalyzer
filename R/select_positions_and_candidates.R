#' Selection of genomic positions
#'
#' Selection of either known methylated positions or suspected methylated positions from data frame (list of data frames) with C-Scores, Z-Scores, GC content information
#'
#' @param wdw.df.list
#' @param annot: a data frame with general project informations. Must have columns specifying the project, rRNA, sample names, paths.
#' @param select.category: a character string specifying the category of methylated positions. Must be one of "known" (default), "suspected" or "candidate"
#' @param annot.name.col: a numerical value specifying the column index of the sample name in the annotation <annot>. Default is 5.
#' @param annot.group.col: a numerical value specifying the column of the biological condition in the annotation <annot>. Default is 3. Note that only 2 different biological conditions are allowed.
#'
#' @return data frame of known (or suspected) methylated positions for all samples
#' @export
#'
#' @examples select.methylated.positions(wdw.df.list, annot, select.category="suspected")
select.methylated.positions <- function(wdw.df.list=NULL, annot=NULL, select.category="known", annot.name.col=5, annot.group.col=3){

  # check that all columns are present
  if (! select.category %in% c("known", "suspected")) {stop("ERROR. Wrong <select.category> parameter. It has to be either <known> or <suspected>.")}
  if (is.null(annot)) {stop("MISSING parameter. Please specify an annotation file <annot>.")}
  if (is.null(wdw.df.list)) {stop("MISSING PARAMETER. Please specify <wdw.df.list>.")}

  all.selected <- lapply(names(wdw.df.list), function(x){

    df          <- wdw.df.list[[x]]
    sample.name <- x

    if (select.category == "known") {
      meth.col <- which((names(df) == "meth") == TRUE) # column with column name <meth>
    } else {
      meth.col <- which((names(df) == "suspected") == TRUE) # column with column name <suspected>
    }

    try({ # exception handling if none of the positions fullfill the conditions which would throw an error when returning empty data frame

        df.meth <- df[which(df[,meth.col] == TRUE),]
        cond <- annot[which(annot[,annot.name.col] == sample.name), annot.group.col]
        df.meth[, "cond"]   <- cond
        df.meth[, "sample"] <- sample.name
        df.meth[, "locus"]  <- c(1:nrow(df.meth))

      }, silent = TRUE)

      if (dim(df.meth)[[1]] == 0) {
        NULL # if emtpy data frame
      } else {
        df.meth # if at least one line in output
      }
  })

  df.all.samples.only.selected            <- do.call("rbind", all.selected)

  if (is.null(df.all.samples.only.selected)) {stop("EMPTY OUTPUT. It seems the the methylation vector is empty (either for known methylated or suspected sites). Please check input parameter <select.col>)")}

  df.all.samples.only.selected[, "cond"]  <- factor(df.all.samples.only.selected[, "cond"], levels=unique(sort(df.all.samples.only.selected[, "cond"])))

  return(df.all.samples.only.selected)
}
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# filters out already known sites
#------------------------------------------------------------------------------
#' Select new candidates based on minimal tresholds for C-Scores and Z-Scores and within a predefined GC-content window. Already known methylated sites are filtered out, suspected sites are kept.
#'
#' @param wdw.df.list
#' @param annot: a data frame with general project informations. Must have columns specifying the project, rRNA, sample names, paths.
#' @param annot.name.col: a numerical value specifying the column index of the sample name in the annotation <annot>. Default is 5.
#' @param annot.group.col: a numerical value specifying the column index of the biological group in the annotation <annot>. Default is 3.
#' @param GC.min: a numeric value between 0 and 1. Minimal GC-content of the surrounding 12-nt window to keep a genomic position and set GC-filter to PASS. Default is 0.2.
#' @param GC.max: a numeric value between 0 and 1. Maximal GC-content of the surrounding 12-nt window to keep a genomic position and set GC-filter to PASS. Default is 0.8.
#' @param zscore: a numeric value above which genomic sites are kept as potential new methylation candidates. Default is 2.
#' @param cscore: a numeric value above which genomic sites are kept as potential new methylation candidates. Default is 0.8.
#'
#' @return Data frame of genomic positions selected to be new methylation candidates corresponding to the criterias (C- and Z-scores, GC-content). Note that positions are kept if they fulfill the 3 criteria in at least one sample of the cohort.
#' @export
#'
#' @examples select.new.candidates(wdw.df.list, annot)
select.new.candidates <- function(wdw.df.list=NULL, annot=NULL, data.position.col=1, annot.name.col=5, annot.group.col=3, GC.min=0.2, GC.max=0.8, zscore=2, cscore=0.8){

  # check that all objects are present
  if (is.null(annot)) {stop("MISSING parameter. Please specify an annotation file <annot>.")}
  if (is.null(wdw.df.list)) {stop("MISSING PARAMETER. Please specify <wdw.df.list>.")}

  all.selected <- lapply(names(wdw.df.list), function(x){

    df          <- wdw.df.list[[x]]
    sample.name <- x

    # Check columns
    if (! "gc" %in% names(df)) {stop("Missing column in input list <wdw.df.list>: column 'gc' for GC content does not exist!")}
    if (! "dist2medInMad" %in% names(df)) {stop("Missing column in input list <wdw.df.list>: column 'dist2medInMad' for Z-score does not exist!")}
    if (! "ScoreC.Median.net" %in% names(df)) {stop("Missing column in input list <wdw.df.list>: column 'ScoreC.Median.net' for C-score does not exist!")}
    if (! "meth" %in% names(df)) {stop("Missing column in input list <wdw.df.list>: column 'meth' for known sites does not exist!")}
    if (! "suspected" %in% names(df)) {stop("Missing column in input list <wdw.df.list>: column 'meth' for suspected sites does not exist!")}

    meth.col      <- which((names(df) == "meth") == TRUE)
    suspected.col <- which((names(df) == "suspected") == TRUE)
    gc.col        <- which((names(df) == "gc") == TRUE)
    cscore.col    <- which((names(df) == "ScoreC.Median.net") == TRUE)
    zscore.col    <- which((names(df) == "dist2medInMad") == TRUE)

    # Filtering
    try({ # exception handling if none of the positions fullfill the conditions which would throw an error when returning empty data frame
      df.filt <- df[which(df[,cscore.col] > cscore),] # filter for minimal C-Score
      df.filt <- df.filt[which(df.filt[,gc.col] > GC.min & df.filt[,gc.col] < GC.max),] # filter based on GC content
      df.filt <- df.filt[which(df.filt[,zscore.col] > zscore),] # filter for minimal z-Score
      df.filt <- df.filt[which(df.filt[,meth.col] != TRUE),] # filter already known sites

      # add columns for condition/sample name/locus
      cond                <- annot[which(annot[,annot.name.col] == sample.name), annot.group.col]
      df.filt[, "cond"]   <- as.factor(cond)
      df.filt[, "sample"] <- sample.name
      #df.filt[, "locus"] <- c(1:nrow(df.filt))

    }, silent = TRUE)

    if (dim(df.filt)[[1]] == 0) {
      NULL # if emtpy data frame
    } else {
      df.filt # if at least one candidate was found
    }

  })

  df.all.samples.new.sites  <- do.call("rbind", all.selected)

  if (is.null(df.all.samples.new.sites)) {stop("EMPTY DATA TABLE. No candidate methylation site selected according to the choosen criteria.")}
  if (! "position" %in% names(df.all.samples.new.sites)) {stop("Missing column in input list <wdw.df.list>: column 'meth' for suspected sites does not exist!")}

  candidates.all.pos        <- as.data.frame(table(df.all.samples.new.sites[, "position"]))

  all.selected.completed <- lapply(names(wdw.df.list), function(x){

    df                      <- wdw.df.list[[x]]
    sample.name             <- x

    candidates              <- df[which(df[,data.position.col] %in% candidates.all.pos[,1]),]
    # add columns for condition/sample name/ locus
    cond                    <- annot[which(annot[,annot.name.col] == sample.name), annot.group.col]
    candidates[, "cond"]    <- as.factor(cond)
    candidates[, "sample"]  <- sample.name
    candidates[, "locus"]   <- c(1:nrow(candidates))
    candidates
  })

  df.all.selected.completed <- do.call("rbind", all.selected.completed)

  return(df.all.selected.completed)

}
#------------------------------------------------------------------------------


