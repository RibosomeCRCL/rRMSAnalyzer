

#-----------------------------------------------------------------------------------------------------------
# Prepares the output data (Cscore and Zscore) for Excel where lignes are positions and columns are samples
# parameters:
# <df.all.samples> : data frame, output of select.methylated.positions()
# <annot> : annotation
# RETURN: returns data frame, 1st column = position and for each sample in df.all.samples 2 columns corresponding to the CScore and to the ZScore
#-----------------------------------------------------------------------------------------------------------
#' Prepares the output data (Cscore and Zscore) for Excel
#'
#' @param df.all.samples
#' @param annot: a data frame with general project informations. Must have columns specifying the project, rRNA, sample names, paths.
#' @param annot.name.col: a numerical value specifying the column of the sample name in the annotation <annot>. Default is 5.
#'
#' @return Returns data frame with 1st column being the position, and two columns for each sample in <df.all.samples> corresponding to the C- and Z-Scores
#' @export
#'
#' @examples prepare.output.data.Scores(df.all.samples, annot)
prepare.output.data.Scores <- function(df.all.samples=NULL, annot=NULL, annot.name.col=5){

  # Check input parameters
  if(is.null(annot)) stop("ERROR: Sample annotation <annot> is missing!")
  #if(is.null(annot.name.col)) stop("ERROR: Column number for <name> column of sample annotation <annot> is missing!")
  if(is.null(df.all.samples)) stop("ERROR: Data is missing! Should be an output of the function: select.methylated.positions().")

  if (sum(c("position", "dist2medInMad", "sample", "ScoreC.Median.net") %in% names(df.all.samples)) != 4) {
    stop("ERROR: Sample data frame column names do not match. At least one column is missing! <df.all.samples> has to be the output of function
            select.methylated.positions() with columns <position>, <dist2medInMad>, <sample>, <ScoreC.Median.net>")
  }


  results.samples <- lapply(annot[, annot.name.col], function(x){

    df.res            <- df.all.samples[which(df.all.samples[, "sample"] == x),]
    col.pos           <- match(c("ScoreC.Median.net", "Z.Score.col.name"), names(df.res))
    df.scores         <- subset(df.res, select=col.pos)
    names(df.scores)  <- c(paste(x, ".C-Score", sep=""), paste(x, ".Z-Score", sep=""))
    df.scores
  })

  Position.Count  <- unique(df.all.samples[, pos.col.name])
  results.excel   <- do.call("cbind", results.samples)
  results.excel   <- cbind(Position.Count, results.excel)
  return(results.excel)
  }
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# Prepares the general informations (annotation) of the output for Excel export, will be completed from the corresponding reference anntation data set
# RETURN: returns data frame with "general" informations: position, gc-content, nomenclatur, NR_Numbering etc
#------------------------------------------------------------------------------
#' Extracts literature information from internal data sets for known/suspected methylation sites
#'
#' Prepares the general informations (annotation) of the output for Excel export, will be completed from the corresponding reference anntation data set parameters:
#'
#' @param df.all.samples: a data frame, output of select.methylated.positions()
#' @param annot: a data frame with general project informations. Must have columns specifying the project, rRNA, sample names, paths.
#' @param select.category: a character string specifying the category of methylated positions. Must be one of "known" (default), "suspected" or "candidate"
#' @param species: a character string specifying the species. Default is "human". In this version of the R-packge, only "human" is supported.
#' @param annot.project.col: a numerical value specifying the column of the project in the annotation <annot>. Default is 1.
#' @param annot.rna.col: a numerical value specifying the column of the rRNA in the annotation <annot>. Default is 4.
#' @param annot.name.col: a numerical value specifying the column of the sample name in the annotation <annot>. Default is 5.
#' @param GC.min: a numeric value between 0 and 1. Minimal GC-content of the surrounding 12-nt window to keep a genomic position and set GC-filter to PASS. Default is 0.2.
#' @param GC.max: a numeric value between 0 and 1. Maximal GC-content of the surrounding 12-nt window to keep a genomic position and set GC-filter to PASS. Default is 0.8.
#' @param verbose:
#'
#' @return Data frame of methylation sites, completed with information from the literature
#' @export
#'
#' @examples general.infos.for.excel.output(df.all.samples, annot, select.category="known", species="human")
general.infos.for.excel.output <- function(df.all.samples=NULL, annot=NULL, select.category="known", species="human", annot.project.col=1, annot.rna.col=4, annot.name.col=5, GC.min=0.2, GC.max=0.8, verbose=TRUE) {

  # Check input parameters
  if(is.null(df.all.samples)) stop("ERROR: Data is missing! Should be output of function: select.methylated.positions().")
  if(is.null(annot)) stop("ERROR: Sample annotation <annot> is missing!")

  if (GC.min > 1 || GC.max > 1 || GC.min < 0 || GC.max < 0) stop("ERROR: GC.min and GC.max must be between 0 and 1!")

  new.candidates <- FALSE

  if (select.category == "known") {
    if (tolower(species) == "human") {
      reference.annot <- human.methylated
      {if (verbose) print("Add information from the literature for KNOWN methylation sites for HUMAN data.")}
    } else if (tolower(species) == "mouse") {
      reference.annot <- mouse.methylated
      {if (verbose) print("Add information from the literature for KNOWN methylation sites for MOUSE data.")}
    } else { stop("ERROR in general.infos.for.excel.output(): Unknown species. Pipeline is only working for human and mouse. Other organisms will be added soon.")}
  } else if (select.category == "suspected") {
    if (tolower(species) == "human") {
       reference.annot <- human.suspected
       {if (verbose) print("Add  information from the literature for SUSPECTED methylation sites for HUMAN data.")}
    } else if (tolower(species) == "mouse") {
        reference.annot <- mouse.suspected
        {if (verbose) print("Add  information from the literature for SUSPECTED methylation sites for MOUSE data.")}
    } else { stop("ERROR in general.infos.for.excel.output(): Unknown species. Pipeline is only working for human and mouse. Other organisms will be added soon.")}
  } else if (select.category == "candidate") {
       new.candidates <- TRUE
       if (verbose) print("No reference annotation provided. Function applied to new methylation candidates.")
  } else { stop("ERROR in general.infos.for.excel.output(): <select.category> parameter has to be one of: known, suspected, candidate!")}


  if(! "sample" %in% names(df.all.samples)) stop("ERROR in general.infos.for.excel.output(): input data frame <df.all.samples> should have a column <sample>")

  if (sum(unique(df.all.samples[, "sample"]) %in% annot[, annot.name.col]) != length(annot[, annot.name.col])) {
    stop("ERROR in general.infos.for.excel.output(): Sample names of input data frame <df.all.samples> do not match those of the sample annotation <annot>.")
  }

  first.sample  <- annot[,annot.name.col][1]

  # take general infos from first.df
  first.df <- df.all.samples[which(df.all.samples[, "sample"] == first.sample),]

  # Check column names
  if (sum(c("position", "type", "identity", "counts", "gc") %in% names(first.df)) != 5) {
    stop("ERROR in general.infos.for.excel.output(): Sample data frame column names do not match. At least one column is missing:
         position, type, identity, counts, gc")
  }

  col.pos.general.infos1      <- match(c("position", "type"), names(first.df))
  general.info1               <- subset(first.df, select=col.pos.general.infos1)

  col.pos.general.infos2      <- match(c("position", "counts", "gc"), names(first.df))
  general.info2               <- subset(first.df, select=col.pos.general.infos2)
  general.info2[, "F.gc_QC"]  <- ifelse(general.info2[, "gc"] < GC.min, "FAIL", ifelse(general.info2[, "gc"] > GC.max, "FAIL", "PASS"))

  col.pos.general.infos3      <- match(c("position", "identity"), names(first.df))
  general.info3               <- subset(first.df, select=col.pos.general.infos3)

  if (new.candidates){ # if function is used for new candidate sites

    if(!((dim(general.info1)[1] == dim(general.info2)[1]) && (dim(general.info2)[1] == dim(general.info3)[1]))) {
      stop("INTERNAL FUNCTION ERROR in general.infos.for.excel.output(): the three general.info[x] data frames do not have same row number.")
    }

    all.infos             <- merge(general.info1, general.info2, by.x = "position", by.y = "position", all.x=T)
    all.infos             <- merge(all.infos, general.info3, by.x = "position", by.y = "position", all.x=T)

    # change column names and orders
    names(all.infos)[1]   <- "Position.Count"
    names(all.infos)[4]   <- "GC.Frequency"

    final.col             <- match(c("type", "Position.Count", "counts", "GC.Frequency", "F.gc_QC", "identity"), names(all.infos))
    all.infos             <- subset(all.infos, select=final.col)
    all.infos.project     <- cbind(project=unique(annot[, annot.project.col]), all.infos)


  } else {

    #rRNA.type            <- unique(annot[, annot.rna.col]) und dann: ref.for.rRNA <- reference.annot[which(reference.annot[, "rRNA"] == r.RNA.type]),]
    ref.for.rRNA          <- reference.annot[which(reference.annot[, "rRNA"] == unique(annot[, annot.rna.col])), ] # select reference annot corresponding to considered rRNA type


    #Check that colnames match
    if (sum(c("Position", "Nomenclature", "NR_046235.Numbering", "SNORD", "Mode.of.coding", "SNORD.host.gene") %in% names(ref.for.rRNA)) != 6) {
      stop("ERROR in general.infos.for.excel.output(): Reference data frame column names do not match. At least one column is missing:
           Position, Nomenclature, NR_046235.Numbering, SNORD, Mode.of.coding, SNORD.host.gene")
    }

    col.pos.ref                 <- match(c("Position", "Nomenclature", "NR_046235.Numbering"), names(ref.for.rRNA))
    part2.names                 <- subset(ref.for.rRNA, select=col.pos.ref)

    col.pos.Snord               <- match(c("Position","SNORD", "Mode.of.coding", "SNORD.host.gene"), names(ref.for.rRNA))
    part4.snord                 <- subset(ref.for.rRNA, select=col.pos.Snord)

    # join ref.for.rRNA and general info
    if(!((dim(general.info1)[1] == dim(general.info2)[1]) && (dim(general.info2)[1] == dim(general.info3)[1]))) {
      stop("INTERNAL FUNCTION ERROR in general.infos.for.excel.output(): the three general.info[x] data frames do not have same row number.")
    }


    if(dim(ref.for.rRNA)[1] != dim(general.info1)[1]) {
      stop("ERROR in general.infos.for.excel.output(): Reference annotation <reference.annot> and selected data frame <general.info1> do not have the same dimensions.
           Some rows are missing. Please check the reference annotation: <human.methylated> for known methylated sites and <suspected.sites> for
           potentially methylated sites.")
    }

    # add project name of annotation file to output data frame

    all.infos             <- merge(general.info1, part2.names, by.x = "position", by.y = "Position", all.x=T)
    all.infos             <- merge(all.infos, general.info2, by.x = "position", by.y = "position", all.x=T)
    all.infos             <- merge(all.infos, part4.snord, by.x = "position", by.y = "Position", all.x=T)
    all.infos             <- merge(all.infos, general.info3, by.x = "position", by.y = "position", all.x=T)

    # change column names and orders
    names(all.infos)[1]   <- "Position.Count"
    names(all.infos)[6]   <- "GC.Frequency"

    final.col <- match(c("type", "Nomenclature", "NR_046235.Numbering", "Position.Count", "counts", "GC.Frequency", "F.gc_QC", "SNORD",
                         "Mode.of.coding", "SNORD.host.gene", "identity"), names(all.infos))

    all.infos             <- subset(all.infos, select=final.col)
    all.infos.project     <- cbind(project=unique(annot[, annot.project.col]), all.infos)

    }

  return(all.infos.project)

}
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
#' Nonparametric Mann-Whitney U-Test to compare methylation between different biological conditions
#'
#' @param output.final
#' @param annot: a data frame with general project informations. Must have columns specifying the project, rRNA, sample names, paths.
#' @param annot.group.col: a numerical value specifying the column of the biological condition in the annotation <annot>. Default is 3. Note that only 2 different biological conditions are allowed.
#' @param annot.name.col: a numerical value specifying the column of the sample name in the annotation <annot>. Default is 5.
#' @param score.to.test: a character string specifying the score to which the Mann-Whitney test is applied. Must be one of "C-Score", "Z-Score", or "both".
#' @param verbose
#' @param digits: an integer value specifying the number of digits for rounding
#'
#' @return Result of statistical test (Nonparametric Mann-Whitney U-Test) for differential methylation between 2 biological conditions
#' @export
#'
#' @examples mann.whitney.test(output.final, annot)
mann.whitney.test <- function(output.final=NULL, annot=NULL, annot.group.col=3, annot.name.col=5, score.to.test="both", verbose=TRUE, digits=5){

  recognize.score <- tolower(sub("\\.", "", sub("-", "", score.to.test))) # check for wrong spelling
  if (recognize.score == "cscore") {
    score.to.test <- "C-Score"
  } else { if (recognize.score == "zscore") {
              score.to.test <- "Z-Score"
  } else {
                if (recognize.score == "both") {
                  score.to.test <- "both"
                } else {
                    stop("ERROR. Wrong <score.to.test> parameter. Must be either <C-Score>, <Z-Score> or <both>.")
                  }
    }
  }

  if (is.null(output.final)) {stop("MISSING parameter. Please specify the data frame <output.final>.")}
  if (is.null(annot)) {stop("MISSING parameter. Please specify an annotation <annot>.")}

  # if test for both scores
  if (tolower(score.to.test) == "both"){
    if (verbose) print("You have chosen to test both C- and Z-scores. First statistical test: C-Score")
    res.CScore      <- mann.whitney.test(output.final=output.final, annot=annot, annot.group.col=annot.group.col, annot.name.col=annot.name.col,
                                                   score.to.test="C-Score")
    if (verbose) print("Second statistical test: Z-Score")
    res.ZScore      <- mann.whitney.test(output.final=output.final, annot=annot, annot.group.col=annot.group.col, annot.name.col=annot.name.col,
                                                   score.to.test="Z-Score")
    res.ZScore.red  <- res.ZScore[, c((ncol(output.final)+1) : ncol(res.ZScore))]

    res.final       <- cbind(res.CScore, res.ZScore.red)
    res.final       <- round.df(res.final, digits=digits)

    return(res.final)
  }

  else {

    groups <- unique(annot[, annot.group.col])
    if (length(groups) != 2) {stop(paste("ERROR: Please check the column with column name '", names(annot)[annot.group.col], "' of your sample annotation <annot>.
                                         You need exactly two conditions."), sep="")}


    all.scores             <- names(output.final)[grep(score.to.test, names(output.final))]
    d                      <- as.matrix(subset(output.final, select=all.scores))

    sample.names           <- colnames(d)

    sample.names           <- substr(sample.names, 1, nchar(sample.names) - nchar(score.to.test) - 1)
    group                  <- as.factor(annot[which(annot[, annot.name.col] == sample.names), annot.group.col]) # get corresponding group for each sample
    mw.pval                <- apply(d, 1, function(row) wilcox.test(row ~ group)$p.value) # Mann whitney U test
    mw.pval.adj            <- p.adjust(mw.pval, method="fdr") # adjust pvalues
    pvals                  <- cbind(mw.pval, mw.pval.adj)
    colnames(pvals)        <- paste(score.to.test, colnames(pvals), sep=".")

    group.mean             <- t(apply(d, 1, function(row) tapply(row, group, mean))) # C-Score group mean
    colnames(group.mean)   <- paste("MEAN", score.to.test, colnames(group.mean), sep=".")

    log2FC                 <- log2(group.mean[,2]/group.mean[,1])
    difference             <- group.mean[,2] - group.mean[,1]
    fold.changes           <- cbind(log2FC, difference)
    colnames(fold.changes) <- paste(score.to.test, c("log2FC", "diff"), sep=".")

    res                   <- cbind(output.final, pvals, group.mean, fold.changes)

    res                   <- round.df(res, digits=digits)
    return(res)
   }

}
#------------------------------------------------------------------------------



#-----------------------------------------------------------------------------------------------------------
# Prepares the output data (Cscore and Zscore) for Excel where lignes are positions and columns are samples
# parameters:
# <df.all.samples> : data frame, output of select.methylated.positions()
# <annot> : annotation
# <annot.name.col> : column number for name in annotation
# <sample.col.name>, default="sample": column name in data frame df.all.samples corresponding to samples
# <Z.Score.col.name>, default = "dist2medInMad": column name in data frame df.all.samples corresponding to Zscore
# <C.Score.col.name>, default = "ScoreC.Median.net": column name in data frame df.all.samples corresponding to Cscore
# pos.col.name, default = "position": column name in data frame df.all.samples corresponding to position
# RETURN: returns data frame, 1st column = position and for each sample in df.all.samples 2 columns corresponding to the CScore and to the ZScore
#-----------------------------------------------------------------------------------------------------------
#' Prepares the output data (C-score and Z-score) for Excel output
#'
#' @param df.all.samples
#' @param annot: a data frame with general project informations. Must have columns specifying the project, rRNA, sample names, paths.
#' @param annot.name.col: a numerical value specifying the column of the sample name in the annotation <annot>. Default is 5.
#'
#' @return data frame with 1st column = position and for each sample in df.all.samples 2 columns corresponding to the CScore and to the ZScore
#' @export
#'
#' @examples prepare.output.data.Scores(df.all.samples, annot)
prepare.output.data.Scores <- function(df.all.samples=NULL, annot=NULL, annot.name.col=5){

  # Check input parameters
  if(is.null(annot)) stop("ERROR: Sample annotation <annot> is missing!")
  #if(is.null(annot.name.col)) stop("ERROR: Column number for <name> column of sample annotation <annot> is missing!")
  if(is.null(df.all.samples)) stop("ERROR: Data is missing! Should be an output of the function: select.methylated.positions().")

 # if (sum(c(pos.col.name, Z.Score.col.name, sample.col.name, C.Score.col.name) %in% names(df.all.samples)) != 4) {
 #  stop("ERROR: Sample data frame column names do not match. At least one column is missing! Check the following parameters:
 #       <pos.col.name>, <sample.col.name>, <Z.Score.col.name>, <C.Score.col.name>")
 #}

  results.samples <- lapply(annot[, annot.name.col], function(x){

    df.res  <- df.all.samples[which(df.all.samples[, "sample"] == x),]
    col.pos <- match(c("ScoreC.Median.net", "dist2medInMad"), names(df.res))
    df.scores <- subset(df.res, select=col.pos)
    names(df.scores) <- c(paste(x, ".C-Score", sep=""), paste(x, ".Z-Score", sep=""))
    df.scores
  }
  )

  Position.Count <- unique(df.all.samples[, "position"])
  results.excel  <- do.call("cbind", results.samples)
  results.excel  <- cbind(Position.Count, results.excel)
  return(results.excel)
}
#------------------------------------------------------------------------------


#' creates output table with scores and information from the literature, if available
#'
#' @param df.all.samples
#' @param annot: a data frame with general project informations. Must have columns specifying the project, rRNA, sample names, paths.
#' @param select.category: a character string specifying the category of methylated positions. Must be one of "known" (default), "suspected" or "candidate"
#' @param species: a character string specifying the species. Default is "human". In this version of the R-packge, only "human" is supported.
#' @param annot.project.col: a numerical value specifying the column of the project in the annotation <annot>. Default is 1.
#' @param annot.rna.col: a numerical value specifying the column of the rRNA in the annotation <annot>. Default is 4.
#' @param annot.name.col: a numerical value specifying the column of the sample name in the annotation <annot>. Default is 5.
#' @param GC.min: a numeric value between 0 and 1. Minimal GC-content of the surrounding 12-nt window to keep a genomic position and set GC-filter to PASS. Default is 0.2.
#' @param GC.max: a numeric value between 0 and 1. Maximal GC-content of the surrounding 12-nt window to keep a genomic position and set GC-filter to PASS. Default is 0.8.
#' @param verbose
#'
#' @return final output data talbe which can be used to compare between conditions
#' @export
#'
#' @examples produce.summary.output.table(df.all.samples, annot, select.category="known", species="human")
produce.summary.output.table <- function(df.all.samples=NULL, annot=NULL, select.category="known", species="human", annot.project.col=1, annot.rna.col=4, annot.name.col=5, GC.min=0.2, GC.max=0.8, verbose=TRUE){

  # Check input parameters
  if(is.null(annot)) stop("ERROR: Sample annotation <annot> is missing!")
  #if(is.null(annot.name.col)) stop("ERROR: Column number for <name> column of sample annotation <annot> is missing!")
  if(is.null(df.all.samples)) stop("ERROR: Data is missing! Should be an output of the function: select.methylated.positions().")

  #print("output.data.scores")
  output.data.scores <- prepare.output.data.Scores(df.all.samples=df.all.samples,
                                                   annot=annot,
                                                   annot.name.col=annot.name.col)
  #print("output.general.infos")
  output.general.infos <- general.infos.for.excel.output(df.all.samples=df.all.samples,
                                                         annot=annot,
                                                         select.category=select.category,
                                                         species=species,
                                                         annot.project.col=annot.project.col,
                                                         annot.rna.col=annot.rna.col,
                                                         annot.name.col=annot.name.col,
                                                         GC.min=GC.min,
                                                         GC.max=GC.max,
                                                         verbose=verbose)
  #print("join")
  output.final.df <- join(output.general.infos, output.data.scores, by="Position.Count")

  return(output.final.df)
}
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
round.df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))

  df[,nums] <- round(df[,nums], digits = digits)

  (df)
}
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
#' Export results to Excel
#'
#' @param results
#' @param output.folder
#' @param select.category
#' @param output.file.prefix
#' @param digits: an interger indicating the number of digits used for rounding the numerical values of the results table
#'
#' @return Exports object to output folder
#' @export
#'
#' @examples export.results.XLS(res, select.category="Methylated")
export.results.XLS <- function(results=NULL, output.folder="", select.category="", output.file.prefix="", digits=5){

  # Check input parameters
  if(is.null(results)) stop("ERROR: The object to be written is missing!")
  #if(is.null(annot.name.col)) stop("ERROR: Column number for <name> column of sample annotation <annot> is missing!")

  if (output.folder == "") {
    output.folder <- getwd()
    print("No output folder provided. Data are written in working directory.")
  }

  # schauen ob liste ist und wieviele elemente
  if (inherits(results, "list")) { # if only 1 data frame
    len <- length(results)
    res.table <- list()
    for (i in 1:len){
      res             <- round.df(results[[i]], digits=digits)
      res.table[[i]]  <- res
    }
    names(res.table) <- names(results)
  } else {
    res.table <- round.df(results, digits=digits)
  }


  # write Excel output
  WriteXLS("res.table", ExcelFileName=paste(output.folder, "/", output.file.prefix, "Scores_", select.category, "_", Sys.Date(), ".xls", sep=""),
           row.names=F)
  print(paste("Results are written to: ", output.folder, "/", output.file.prefix, "Scores_", select.category, "_", Sys.Date(), ".xls", sep=""))

}

