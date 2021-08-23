
#--------------------------------------------------------------------------------
# compute C-Score based on median
#--------------------------------------------------------------------------------

#' Computation of C- and Z-scores for genomic positions within a predefined window
#'
#' @param ds
#' @param known.meth.vector
#' @param suspected.meth.vector
#' @param flanking
#' @param data.counts.col
#' @param data.position.col
#'
#' @return data frame with C-scores and Z-scores for each genomic postion
#' @export compute.window.median.mad
#'
#' @examples compute.window.median.mad(ds, known.meth.vector, suspected.meth.vector)
compute.window.median.mad <- function(ds=NULL, known.meth.vector=NULL, suspected.meth.vector=NULL, flanking=6, data.counts.col=4, data.position.col=1) {

  # check that all parameters exist
  if (is.null(ds)) {stop("MISSING parameter. Please specify a data frame <ds>.")}
  if (is.null(known.meth.vector)) {stop("MISSING PARAMETER in compute.window.median.mad(). Please specify <known.meth.vector>.")}
  if (is.null(suspected.meth.vector)) {stop("MISSING PARAMETER in compute.window.median.mad(). Please specify <suspected.meth.vector>.")}

  if (! is.null(known.meth.vector)) {
    ds[, "meth"] <- ds[, data.position.col] %in% known.meth.vector
  }
  if (! is.null(suspected.meth.vector)) {
    ds[, "suspected"] <- ds[, data.position.col] %in% suspected.meth.vector
  }

  ds[, "median"] <- NA
  ds[, "mad"] <- NA

  # compute scores
  for (i in (flanking + 1) : (dim(ds)[1] - flanking)){

    ds[i, "median"] <- median(ds[c((i-flanking) : (i-1), (i+1) : (i + flanking)), data.counts.col]) # median
    ds[i, "mad"] <- mad(ds[c((i-flanking) : (i-1), (i+1) : (i + flanking)), data.counts.col]) # mad

  }
  ds[, "dist2medInMad"] <- (ds[,"median"] - ds[, data.counts.col])/ds[,"mad"]
  scoreC.Median.raw <- 1 - ds[, data.counts.col]/ds[, "median"]
  ds[, "ScoreC.Median.net"] <- ifelse(scoreC.Median.raw < 0, 0, scoreC.Median.raw)

  return(ds)
}
#------------------------------------------



#------------------------------------------
# Compute Cscore and Zscore
#------------------------------------------
#' Computes C- and Z-scores of all samples of a given study
#'
#' @param all.samples
#' @param annot: a data frame with general project informations. Must have columns specifying the project, rRNA, sample names, paths.
#' @param annot.rna.col: a numerical value specifying the column of the rRNA in the annotation <annot>. Default is 4.
#' @param species: a character string specifying the species. Default is "human". In this version of the R-packge, only "human" is supported.
#' @param window.size
#' @param data.position.col
#' @param data.counts.col
#' @param gc.probas.pos.col
#' @param gc.probas.GC.col
#'
#' @return List of data frames. Each data frame is a sample (count table) for which the C- and Z-scores have been computed
#' @export compute.cscore.zscore
#'
#' @examples compute.cscore.zscore(all.samples, annot)
compute.cscore.zscore <- function(all.samples=NULL, annot=NULL, annot.rna.col=4, species = "human", verbose=FALSE, window.size=6, data.position.col=1, data.counts.col=4, gc.probas.pos.col=1, gc.probas.GC.col=7){

  if (is.null(all.samples)) {stop("MISSING parameter. Please specify <all.samples>.")}
  if (is.null(annot)) {stop("MISSING parameter. Please specify <all.samples>.")}
  if (is.null(window.size)) {stop("MISSING PARAMETER. Please specify the window size for the score computation, <window.size>.")}
  if (is.null(data.position.col)) {stop("MISSING PARAMETER. Please specify <data.position.col>.")}
  #if ((gc.probas.pos.col)): --> position colum in  <wdw.probas>, internal data frame, used in <load.gc.content.for.human.ribosome()>
  #if ((gc.probas.GC.col)): -->  GC column in <wdw.probas>, internal data frame, used in <load.gc.content.for.human.ribosome()>


  if (tolower(species) == "human") {
    gc.rna <- load.gc.content.for.human.ribosome(annot=annot, annot.rna.col=annot.rna.col)
  } else if (tolower(species) == "mouse"){
    gc.rna <- load.gc.content.for.mouse.ribosome(annot=annot, annot.rna.col=annot.rna.col)
  } else {stop("ERROR. Unknown species. Pipeline is only working for human and mouse. Other organisms will be added soon.")}

  pos.known.meth      <- gc.rna[[2]]
  pos.suspected.meth  <- gc.rna[[3]]
  gc.probas           <- gc.rna[[1]]

  wdw.df.list <- lapply(all.samples, function(df){

    gc <- gc.probas[match(df[, data.position.col], gc.probas[, gc.probas.pos.col]), gc.probas.GC.col ]
    df <- cbind(df, gc)

    df2 <- compute.window.median.mad(ds                    = df,
                                     flanking              = window.size,
                                     data.counts.col       = data.counts.col,
                                     data.position.col     = data.position.col,
                                     known.meth.vector     = pos.known.meth,
                                     suspected.meth.vector = pos.suspected.meth)
    df2
  })
  return(wdw.df.list)
}
#------------------------------------------

#' Calculate both C-score and Z-score for one RNA of one sample.
#'
#' @param ds 
#' @param flanking 
#' @param data.counts.col 
#' @param data.position.col 
#'
#' @return
#' @export
#'
#' @examples
calculate_score_by_RNA <- function(ds, flanking=6, data_counts_col=4, data_position_col=1) {
  
  # check that all parameters exist
  if (is.null(ds)) {stop("MISSING parameter. Please specify a data frame <ds>.")}
  ds[, "mean"] <- NA
  ds[, "median"] <- NA
  ds[, "mad"] <- NA
  
  # compute scores
  for (i in (flanking + 1) : (dim(ds)[1] - flanking)){
    #TODO : selection mean/median
    ds[i, "mean"] <- mean(ds[c((i-flanking) : (i-1), (i+1) : (i + flanking)), data_counts_col]) # mean
    ds[i, "median"] <- median(ds[c((i-flanking) : (i-1), (i+1) : (i + flanking)), data_counts_col]) # median
    ds[i, "mad"] <- mad(ds[c((i-flanking) : (i-1), (i+1) : (i + flanking)), data_counts_col]) # mad
    
  }
  ds[, "dist2medInMad"] <- (ds[,"median"] - ds[, data_counts_col])/ds[,"mad"]
  scorec_median_raw <- 1 - ds[, data_counts_col]/ds[, "median"]
  ds[, "ScoreC.Median.net"] <- ifelse(scorec_median_raw < 0, 0, scorec_median_raw)
  
  return(ds)
}


#' Title
#'
#' @param sample_df 
#' @param data_rna_col 
#' @param flanking 
#' @param data_counts_col 
#' @param data_position_col 
#'
#' @return
#' @export
#'
#' @examples
calculate_score_by_sample <- function(sample_df=NULL,
                                      data_rna_col = 1,
                                      flanking=6, 
                                      data_counts_col=3,
                                      data_position_col=1) {
  
  sample_df[,data_rna_col] <- as.factor(sample_df[,data_rna_col])
  RNA_counts_list <- split(sample_df, sample_df[,data_rna_col])
  sample_score <- lapply(RNA_counts_list, calculate_score_by_RNA, flanking = flanking, data_counts_col = data_counts_col, data_position_col = data_position_col)
  
  return(dplyr::bind_rows(sample_score))
  
}

#' Calculate score for a whole sample cohort
#'
#' @param dt 
#' @param flanking 
#' @param data_counts_col 
#' @param data_position_col 
#'
#' @return
#' @export
#'
#' @examples
calculate_score <- function(dt=NULL, flanking=6,
                            data_rna_col = 1,
                            data_counts_col=4,
                            use_multithreads = F) {
  
  if(class(dt) == "RiboClass") {
    dt = dt["raw_counts"] #we only need the counts to calculate the score
  }
   # Experimental : Multithreading is 3x faster than single-thread
   # TODO : implement 
   if(use_multithreads) samples_czscore <- BiocParallel::bplapply(dt[["raw_counts"]], calculate_score_by_sample)
   else samples_czscore <- lapply(dt[["raw_counts"]], calculate_score_by_sample, data_counts_col = data_counts_col, data_rna_col = data_rna_col)
  
  return(samples_czscore)
  
  
}
  
