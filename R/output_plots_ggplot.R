#' C-score plot along the genomic positions
#'
#'
#' @param df: a data frame to be plotted. Data frame must have the following columns: position, type, identity, ScoreC.Median.net, cond, sample, locus.
#' @param save.as.png: a logical indicating if the graphical output should be saved as .png or not. Default is TRUE.
#' @param display.plot: a logical indicating if the graphical output should be displayed at the default screen device. Default is TRUE.
#' @param output.folder: a character string indicating a path to an existing folder where plots are saved. Default is the working directory.
#' @param plot.analysis.type: a character string used to construct the output file name. Suggestions: "Methylated", "Suspected", "Candidates". Default is "UNDEFINED".
#' @param plot.new.candidates: a logical (default is FALSE) indicating if plot is created for data from known/suspected methylation sites or to plot new candidate sites in order to define the label of facets. The identity column is used for known/suspected sites. The rRNA_position for new candidate sites.
#' @param plot.prefix: a character string used to construct the output file names. Note that the output file will start with this string.
#' @param verbose
#'
#' @return C-score profiles along the genomic positions of the data frame. Biological conditions appear as different panels. Samples are diplayed by different colors.
#' @export plot.cscores.by.site
#'
#' @examples plot.cscores.by.site(df)
plot.cscores.by.site <- function(df=NULL, save.as.png=TRUE, display.plot=TRUE, output.folder="", plot.analysis.type="UNDEFINED", plot.new.candidates=FALSE, plot.prefix="", verbose=TRUE){

  if(is.null(df)) stop("ERROR: Missing input data set <df.all.samples.selected.pos>")

  if (output.folder == "") {
    output.folder <- getwd()

    if (verbose) print("No output folder provided. Images are saved in working directory.")
  }

  #------------------------------------------
  # Parametrage for pipeline
  #------------------------------------------
  res.col.names <- names(df)

  # check that all columns needed are present
  if (! "cond" %in% res.col.names) {stop("Missing column in input data frame <df>: column 'cond' does not exist!")}
  if (! "sample" %in% res.col.names) {stop("Missing column in input data frame <df>: column 'sample' does not exist!")}
  if (! "locus" %in% res.col.names) {stop("Missing column in input data frame <df>: column 'locus' does not exist!")}
  if (! "identity" %in% res.col.names) {stop("Missing column in input data frame <df>: column 'identity' does not exist!")}
  if (! "ScoreC.Median.net" %in% res.col.names) {stop("Missing column in input data frame <df>: column 'ScoreC.Median.net' does not exist!")}
  if (! "position" %in% res.col.names) {stop("Missing column in input data frame <df>: column 'position' does not exist!")}
  if (! "type" %in% res.col.names) {stop("Missing column in input data frame <df>: column 'type' does not exist!")}

  #res.cond.col   <- which(res.col.names == "cond")
  #res.sample.col <- which(res.col.names == "sample")
  #res.locus.col  <- which(res.col.names == "locus")
  res.id.col     <- which(res.col.names == "identity")
  #res.ScoreC.col <- which(res.col.names == "ScoreC.Median.net")
  res.pos.col    <- which(res.col.names == "position")
  res.type.col   <- which(res.col.names == "type")


  #------------------------------------------
  # Plots
  #------------------------------------------

  if (plot.new.candidates) {
    if (verbose) print("C-score plot for NEW CANDIDATES")
    df[, res.id.col] <- paste(df[, res.type.col], df[, res.pos.col], sep="_")
    column.for.plot <- res.id.col # if new candidates and thus identity column = NA for almost all sites; vorher: res.pos.col
  } else {
    if (verbose) print("C-score plot for KNOWN or SUSPECTED sites")
    column.for.plot <- res.id.col
  }

  # change xlabels
  breaks <- c(1:length(unique(df[,column.for.plot])))
  xlabels <- unique(df[,column.for.plot])

  # define colors
  all.combinations <- table(df$sample, df$cond)
  nb.cond <- dim(all.combinations)[2]
  nb.samp <- max(table(all.combinations[,1]))

  cscores.gg <- ggplot(data=df, aes(x=factor(locus), y=ScoreC.Median.net, group=sample, colour=sample)) +
                  scale_x_discrete(breaks=breaks, labels=xlabels, name="rRNA methylation sites") +
                  geom_point() +
                  geom_line() +
                  facet_wrap(~cond, ncol=1) +
                  theme_bw() +
                  theme(axis.text.x = element_text(angle=90)) +
                  xlab("rRNA methylation sites") +
                  ylab("C-Score (Median)") +
                  ylim(0,1) +
                  #scale_colour_manual(values=rep(viridis(nb.samp), nb.cond)) # for viridis color scale
                  scale_colour_manual(values=rep(hue_pal()(nb.samp), nb.cond))

  if (save.as.png) {
    if (verbose) print(paste("xyplot: ", output.folder, "/", plot.prefix, "QC_3_",  plot.analysis.type, "_C-Score_profiles_", Sys.Date(), ".png", sep=""))
    png(file=paste(output.folder, "/", plot.prefix, "QC_3_",  plot.analysis.type, "_C-Score_profiles_", Sys.Date(), ".png", sep=""), height = 600, width = 1000)
      print(cscores.gg)
    dev.off()
  }
  if (display.plot){
    print(cscores.gg)
  }

}
#----------------------------------------------------------------------------------






#----------------------------------------------------------------------------------
#' Boxplot of C-scores across biological conditions
#'
#' @param df: a data frame to be plotted. Data frame must have the following columns: position, type, identity, ScoreC.Median.net, cond, sample, locus.
#' @param save.as.png: a logical indicating if the graphical output should be saved as .png or not. Default is TRUE.
#' @param display.plot: a logical indicating if the graphical output should be displayed at the default screen device. Default is TRUE.
#' @param output.folder: a character string indicating a path to an existing folder where plots are saved. Default is the working directory.
#' @param plot.analysis.type: a character string used to construct the output file name. Suggestions: "Methylated", "Suspected", "Candidates". Default is "UNDEFINED".
#' @param plot.new.candidates: a logical (default is FALSE) indicating if plot is created for data from known/suspected methylation sites or to plot new candidate sites in order to define the label of facets. The identity column is used for known/suspected sites. The rRNA_position for new candidate sites.
#' @param plot.prefix: a character string used to construct the output file names. Note that the output file will start with this string.
#' @param verbose
#'
#' @return Boxplot of C-scores by different biological groups. Each genomic positions ons appear as different panels. Samples are diplayed by different colors.
#' @export plot.cscores.boxplot
#'
#' @examples plot.cscores.boxplot(df)
plot.cscores.boxplot <- function(df=NULL, save.as.png=TRUE, display.plot=TRUE, output.folder="", plot.analysis.type="UNDEFINED", plot.new.candidates=FALSE, plot.prefix="", verbose=TRUE){

  if(is.null(df)) stop("ERROR: Missing input data set <df.all.samples.selected.pos>")

  if (output.folder == "") {
    output.folder <- getwd()
    if (verbose) print("No output folder provided. Images are saved in working directory.")
  }

  #------------------------------------------
  # Parametrage for pipeline
  #------------------------------------------
  res.col.names <- names(df)

  # check that all columns needed are present
  if (! "cond" %in% res.col.names) {stop("Missing column in input data frame <df>: column 'cond' does not exist!")}
  if (! "sample" %in% res.col.names) {stop("Missing column in input data frame <df>: column 'sample' does not exist!")}
  if (! "locus" %in% res.col.names) {stop("Missing column in input data frame <df>: column 'locus' does not exist!")}
  if (! "identity" %in% res.col.names) {stop("Missing column in input data frame <df>: column 'identity' does not exist!")}
  if (! "ScoreC.Median.net" %in% res.col.names) {stop("Missing column in input data frame <df>: column 'ScoreC.Median.net' does not exist!")}
  if (! "position" %in% res.col.names) {stop("Missing column in input data frame <df>: column 'position' does not exist!")}
  if (! "type" %in% res.col.names) {stop("Missing column in input data frame <df>: column 'type' does not exist!")}

  res.cond.col   <- which(res.col.names == "cond") # used
  res.sample.col <- which(res.col.names == "sample")
  res.locus.col  <- which(res.col.names == "locus")
  res.id.col     <- which(res.col.names == "identity")
  res.ScoreC.col <- which(res.col.names == "ScoreC.Median.net")
  res.pos.col    <- which(res.col.names == "position")
  res.type.col   <- which(res.col.names == "type")


  if (plot.new.candidates) {
    if (verbose) print("C-score plot for NEW CANDIDATES")
    df[, res.id.col] <- paste(df[, res.type.col], df[, res.pos.col], sep="_")
    column.for.plot <- res.id.col # if new candidates and thus identity column = NA for almost all sites; vorher: res.pos.col
  } else {
    if (verbose) print("C-score plot for KNOWN or SUSPECTED sites")
    column.for.plot <- "identity"
  }

  ncol.gg <- ceiling(sqrt(length(unique(df[,res.locus.col])))) # number of columns to get squared plot

  # to get panel names in order of the data frame
  df$id.new <- factor(df[, res.id.col], levels = unique(df[, res.id.col]), labels=unique(df[, res.id.col]))

  # compute pvalues: Mann-Whitney-U test with correction for multiple test
  pvals1 <- lapply(unique(df[, res.id.col]), function(x) df[ which(df[, res.id.col] == x),])
  names(pvals1) <- unique(df[, res.id.col])
  pvals2 <- lapply(pvals1, function(x){
              mu.test <- wilcox.test(x[, res.ScoreC.col] ~ x[,res.cond.col])$p.value
              mu.test
  })
  pvals <- data.frame(id.new = names(pvals2), pval= paste("p=", round(p.adjust(unlist(pvals2), method="fdr"), digits=3), sep=""))

  #boxplot
  boxplot.gg <- ggplot(data=df, aes(x=cond, y=ScoreC.Median.net, group=cond, colour=cond)) +
                  geom_boxplot() +
                  facet_wrap(~id.new, ncol=ncol.gg) +
                  theme_bw() +
                  xlab("") +
                  ylab("C-Score (Median)") +
                  theme(axis.text.x = element_text(angle=90), legend.position="bottom", legend.title=element_blank()) +
                  ylim(0,1) +
                  #stat_compare_means(label.y=0.1, label="p.format", hide.ns=TRUE)
                  geom_text(data=pvals, aes(x=-Inf, y=-Inf, label=pval), hjust=-0.1, vjust=-1, inherit.aes = FALSE, size=3)

                  #compare_means(ScoreC.Median.net ~ cond, df, group.by="id.new", p.adjust.method = "fdr")

  if (save.as.png) {
    if (verbose) print(paste("boxplot: ", output.folder, "/", plot.prefix, "Results_2_", plot.analysis.type, "_C-Score_boxplot_", Sys.Date(), ".png", sep=""))
    png(file=paste(output.folder, "/", plot.prefix, "Results_2_", plot.analysis.type, "_C-Score_boxplot_", Sys.Date(), ".png", sep=""), height = 600, width = 1000)
      print(boxplot.gg)
    dev.off()
  }
  if (display.plot){
    print(boxplot.gg)
  }

}
#----------------------------------------------------------------------------------



#----------------------------------------------------------------------------------
#' C-scores by biological conditions across all samples
#'
#' Graphical representation of the C-Scores. Each panel corresponds to one sample. Represented are either known, suspected or new methylation sites
#'
#' @param df: a data frame to be plotted. Data frame must have the following columns: position, type, identity, ScoreC.Median.net, cond, sample, locus.
#' @param save.as.png: a logical indicating if the graphical output should be saved as .png or not. Default is TRUE.
#' @param display.plot: a logical indicating if the graphical output should be displayed at the default screen device. Default is TRUE.
#' @param output.folder: a character string indicating a path to an existing folder where plots are saved. Default is the working directory.
#' @param plot.analysis.type: a character string used to construct the output file name. Suggestions: "Methylated", "Suspected", "Candidates". Default is "UNDEFINED".
#' @param plot.new.candidates: a logical (default is FALSE) indicating if plot is created for data from known/suspected methylation sites or to plot new candidate sites in order to define the label of facets. The identity column is used for known/suspected sites. The rRNA_position for new candidate sites.
#' @param plot.prefix: a character string used to construct the output file names. Note that the output file will start with this string.
#' @param verbose
#'
#' @return plot
#' @export plot.cscores.dotplot
#'
#' @examples plot.cscores.dotplot(df)
plot.cscores.dotplot <- function(df=NULL, save.as.png=TRUE, display.plot=TRUE, output.folder="", plot.analysis.type="UNDEFINED", plot.new.candidates=FALSE, plot.prefix="", verbose=TRUE){

  if(is.null(df)) stop("ERROR: Missing input data set <df.all.samples.selected.pos>")

  if (output.folder == "") {
    output.folder <- getwd()
    if (verbose) print("No output folder provided. Images are saved in working directory.")
  }

  #------------------------------------------
  # Parametrage for pipeline
  #------------------------------------------
  res.col.names <- names(df)

  # check that all columns needed are present
  if (! "cond" %in% res.col.names) {stop("Missing column in input data frame <df>: column 'cond' does not exist!")}
  if (! "sample" %in% res.col.names) {stop("Missing column in input data frame <df>: column 'sample' does not exist!")}
  if (! "locus" %in% res.col.names) {stop("Missing column in input data frame <df>: column 'locus' does not exist!")}
  if (! "identity" %in% res.col.names) {stop("Missing column in input data frame <df>: column 'identity' does not exist!")}
  if (! "ScoreC.Median.net" %in% res.col.names) {stop("Missing column in input data frame <df>: column 'ScoreC.Median.net' does not exist!")}
  if (! "position" %in% res.col.names) {stop("Missing column in input data frame <df>: column 'position' does not exist!")}
  if (! "type" %in% res.col.names) {stop("Missing column in input data frame <df>: column 'type' does not exist!")}

  res.cond.col   <- which(res.col.names == "cond")
  res.sample.col <- which(res.col.names == "sample")
  res.locus.col  <- which(res.col.names == "locus")
  res.id.col     <- which(res.col.names == "identity")
  res.ScoreC.col <- which(res.col.names == "ScoreC.Median.net")
  res.pos.col    <- which(res.col.names == "position")
  res.type.col   <- which(res.col.names == "type")


  if (plot.new.candidates) {
    if (verbose) print("C-score plot for NEW CANDIDATES")
    df[, res.id.col] <- paste(df[, res.type.col], df[, res.pos.col], sep="_")
    column.for.plot <- res.id.col # if new candidates and thus identity column = NA for almost all sites; vorher: res.pos.col
  } else {
    if (verbose) print("C-score plot for KNOWN or SUSPECTED sites")
    column.for.plot <- "identity"
  }

  ncol.gg <- ceiling(sqrt(length(unique(df[,res.locus.col])))) # number of columns to get squared plot

  df$id.new <- factor(df[, res.id.col], levels = unique(df[, res.id.col]), labels=unique(df[, res.id.col]))

  dotplot.gg <- ggplot(data=df, aes(x=cond, y=ScoreC.Median.net, group=cond, colour=cond, label = sample)) +
                  geom_point(aes(fill=cond)) +
                  facet_wrap(~id.new, ncol=ncol.gg) +
                  theme_bw() +
                  xlab("") +
                  ylab("C-Score (Median)") +
                  theme(axis.text.x = element_text(angle=90),legend.title=element_blank()) +
                  ylim(0,1)


  if (save.as.png) {
    if (verbose) print(paste("dotplot: ", output.folder, "/", plot.prefix, "Results_1_", plot.analysis.type, "_C-Score_dotplot_", Sys.Date(), ".png", sep=""))
    png(file=paste(output.folder, "/", plot.prefix, "Results_1_", plot.analysis.type, "_C-Score_dotplot_", Sys.Date(), ".png", sep=""), height = 600, width = 1000)
      print(dotplot.gg)
    dev.off()
  }

  if (display.plot){
    print(dotplot.gg)
  }

}
#------------------------------------------------------------------------------




#----------------------------------------------------------------------------------

#' C-score versus Z-score
#'
#' @param wdw.df.list
#' @param save.as.png: a logical indicating if the graphical output should be saved as .png or not. Default is TRUE.
#' @param display.plot: a logical indicating if the graphical output should be displayed at the default screen device. Default is TRUE.
#' @param output.folder: a character string indicating a path to an existing folder where plots are saved. Default is the working directory.
#' @param plot.prefix: a character string used to construct the output file names. Note that the output file will start with this string.
#' @param xyplot.rows
#' @param zscore.col.name
#' @param cscore.col.name
#' @param meth.col.name
#' @param suspected.col.name
#' @param verbose
#' @param zscore.threshold: a numeric value indicating the Z-score threshold for the plot, which is represented by a dotted line. Default is 2.
#' @param cscore.threshold: a numeric value indicating the C-score threshold for the plot, which is represented by a dotted line. Default is 0.8.
#'
#' @return plot
#' @export plot.cscore.vs.zscore.gg
#'
#' @examples plot.cscore.vs.zscore.gg(wdw.df.list)
plot.cscore.vs.zscore.gg <- function(wdw.df.list=NULL, save.as.png=TRUE, display.plot=TRUE, output.folder="", plot.prefix="", xyplot.rows=2, zscore.col.name="dist2medInMad", cscore.col.name="ScoreC.Median.net", meth.col.name="meth", suspected.col.name="suspected", verbose=TRUE, zscore.threshold=2, cscore.threshold=0.8 ){

  # Check input parameters
  if(is.null(wdw.df.list)) stop("ERROR: Input data <wdw.df.list> is missing!")
  if(is.null(output.folder)) stop("ERROR: Please set the path <output.folder> where the images will be saved.")
  if(is.null(plot.prefix)) stop("ERROR: <plot.prefix> for plot titles is missing.")

  if (output.folder == "" && save.as.png) {
    output.folder <- getwd()
    print("No output folder provided. Images are saved in working directory.")
  }

  xyplot.cols <- ceiling(length(wdw.df.list)/xyplot.rows)

  df.wdw <- do.call("rbind", wdw.df.list)
  s.id <- sapply(rownames(df.wdw), function(x) {
            l <- strsplit(x, "\\.")[[1]][1]
            l
  })

  if (! "dist2medInMad" %in% names(df.wdw)) {stop("Missing column in input data frame <wdw.df.list>: column 'dist2medInMad' does not exist!")}
  if (! "ScoreC.Median.net" %in% names(df.wdw)) {stop("Missing column in input data frame <wdw.df.list>: column 'ScoreC.Median.net' does not exist!")}
  if (! "meth" %in% names(df.wdw)) {stop("Missing column in input data frame <wdw.df.list>: column 'meth' does not exist!")}
  if (! "suspected" %in% names(df.wdw)) {stop("Missing column in input data frame <wdw.df.list>: column 'suspected' does not exist!")}

  df.wdw[, "sampleID"]    <- as.factor(s.id)
  df.wdw[, "meth.status"] <- as.factor(ifelse(df.wdw[, "meth"] == T, "methylated", ifelse(df.wdw[, "suspected"] == T, "suspected", "others")))

  #http://rstudio-pubs-static.s3.amazonaws.com/250110_6dbc9e139a8e461a85ff3d017e8565fd.html -> code for transparency
  zscore.cscore.plot <- ggplot(data=df.wdw,
                            aes(y=ScoreC.Median.net, x=dist2medInMad, color=meth.status, alpha=meth.status)) +
                            facet_wrap(~sampleID, ncol=xyplot.cols) +
                            geom_point() +
                            ylab("C-Score (Median)") +
                            xlab("Z-score") +
                            xlim(0,6) +
                            ylim(0,1) +
                            theme_bw() +
                            geom_vline(xintercept=zscore.threshold, linetype="dotted") +
                            geom_hline(yintercept=cscore.threshold, linetype="dotted") +
                            scale_color_manual(values=c("red", "black", "green")) +
                            scale_alpha_manual(values=c(0.6, 1, 0.9), guide=F)


  if (save.as.png) {
    if (verbose) print(paste(output.folder, "/", plot.prefix, "Zscore_vs_Cscore_", Sys.Date(), ".png", sep=""))
    png(file=paste(output.folder, "/", plot.prefix, "Zscore_vs_Cscore_", Sys.Date(), ".png", sep=""), height = 600, width = 1000)
      print(zscore.cscore.plot)
    dev.off()
  }
  if (display.plot){
    print(zscore.cscore.plot)
  }

}
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
#' Volcano plot of C-Scores
#'
#' @param result.statistical.test
#' @param score.to.test: a string indicating the type of score to use for the volcano plot. Can by <C-Score> or <Z-Score>. Default is <C-Score>.
#' @param save.as.png: a logical indicating if the graphical output should be saved as .png or not. Default is TRUE.
#' @param output.folder: a character string indicating a path to an existing folder where plots are saved. Default is the working directory.
#' @param plot.prefix: a character string used to construct the output file names. Note that the output file will start with this string.
#' @param plot.analysis.type: a character string used to construct the output file name. Suggestions: "Methylated", "Suspected", "Candidates". Default is "UNDEFINED".
#' @param display.plot: a logical indicating if the graphical output should be displayed at the default screen device. Default is TRUE.
#' @param verbose
#' @param xintercept.neg: a numerical value indicating an x-value where to add a dotted vertical line at the volcano plot.
#' @param xintercept.pos: a numerical value indicating an x-value where to add a dotted vertical line at the volcano plot.
#' @param y.intercept: a numerical value indicating an y-value where to add a dotted horizontal line at the volcano plot. All points beyond will be displayed with their labels
#'
#' @return Volcano plot (png) of scores, using fold change and difference on x-axis
#' @export plot.volcano.of.scores
#' @usage plot.volcano.of.scores
#'
#' @examples plot.volcano.of.scores(result.statistical.test, output.folder)
plot.volcano.of.scores <- function(result.statistical.test=NULL, save.as.png=TRUE, display.plot=TRUE, output.folder="", plot.prefix="", plot.analysis.type="UNDEFINED", score.to.test="C-Score", verbose=TRUE, xintercept.neg=-0.2, xintercept.pos=0.2, y.intercept=1.3){

  if(is.null(result.statistical.test)) stop("ERROR: Missing input data set <result.statistical.test>")
  if(is.null(plot.analysis.type)) stop("ERROR: Missing input data parameter <plot.analysis.type>. Should be <Methylated>, <Suspected> or <New_candidates>")

  if (output.folder == "") {
    output.folder <- getwd()
    print("No output folder provided. Images are saved in working directory.")
  }


  mean.columns <- grep("MEAN", names(result.statistical.test))
  pval.column  <- grep("pval.adj", names(result.statistical.test))
  id.column    <- grep("identity", names(result.statistical.test))
  type.column  <- grep("type", names(result.statistical.test))
  pos.column   <- grep("Position.Count", names(result.statistical.test))

  if (length(mean.columns) < 2) stop("ERROR: Please check your input data <result.statistical.test>. There have to be exactly 2 columns with the
                                     group means")
  # For data frames that computed both Z-Score and C-Score --> select columns according to <score.to.test> parameter
  if (length(mean.columns) > 2) {

    # Hier noch irgendwo hin schreiben, welchen Score man testet
    #if (is.null(score.to.test)) stop("ERROR: Your input data has pvalues for both C-Score and Z-Score. Please specify <score.to.test>!")
    if (! score.to.test %in% c("C-Score", "Z-Score")) {stop("ERROR. Wrong <score.to.test> parameter. <score.to.test> has to be either <C-Score>, <Z-Score>.")}

    mean.columns.test <- grep(score.to.test, names(result.statistical.test)[mean.columns])
    mean.columns      <- names(result.statistical.test)[mean.columns][mean.columns.test]

    if (length(mean.columns) < 2) stop(paste("ERROR: Please check your input data <result.statistical.test>. There have to be exactly 2 columns with the
                                             group means for the selected score <score.to.test> ", score.to.use, sep=""))
    pval.column.test <- grep(score.to.test, names(result.statistical.test)[pval.column])
    pval.column      <- names(result.statistical.test)[pval.column][pval.column.test]

  }


  # Fold change for volcano plot
  fc <- log2(result.statistical.test[, mean.columns[2]] / result.statistical.test[, mean.columns[1]])
  fc.x.axis <- sapply(mean.columns, function(x) substr(x, 14, nchar(x))) # note that 14 = length of <MEAN.C-Score.> or <MEAN.Z-Score.>

  df <- data.frame(fc=fc, pval=-log10(result.statistical.test[, pval.column]))

  volcano.plot.fc <- ggplot(data=df) + geom_point(aes(x=fc, y=pval), size=2, shape=20) +
          ylab("-log10(adjusted p-value)") +
          xlab(paste("log2(fold change), fold change =", fc.x.axis[2], "/", fc.x.axis[1], sep="")) +
          theme_bw() +
          geom_vline(xintercept=0.05, linetype="dotted") +
          geom_vline(xintercept=-0.05, linetype="dotted") +
          geom_hline(yintercept=1.3, linetype="dotted")

  # Difference
  df2 <- data.frame(diff =result.statistical.test[, mean.columns[2]] - result.statistical.test[, mean.columns[1]],
                    pval=-log10(result.statistical.test[, pval.column]))

  if( sum(is.na(result.statistical.test[,id.column])) == 0){
    df2[, "labels"] <- result.statistical.test[,id.column]
  } else {
    df2[, "labels"] <- paste(result.statistical.test[, type.column], result.statistical.test[, pos.column], sep="_")
  }
  df2[, "status"]   <- as.factor((ifelse(df2[, "pval"] > y.intercept, 2, 1)))
  
  max.pval <- max(df2$pval)


  volcano.plot.diff <- ggplot(data=df2, aes(x=diff, y=pval, label=labels, color=status)) +
                       geom_point(size=2, shape=20) +
                       # geom_text(data=df2, aes_string(label=ifelse(status == 2, as.character(labels), '')), hjust=-0.1, vjust=1, angle=90, color="red") +
                       # geom_text_repel(data=subset(df2, status > 1), aes(diff, pval, label=labels, color="red"), color="red") +
                       geom_text_repel(data=subset(df2, status == 2), aes(diff, pval, label=labels, color="red"), color="red") +
                       ylab("-log10(adjusted p-value)") +
                       xlab(paste("Difference (", fc.x.axis[2], "-", fc.x.axis[1], ")", sep="")) +
                       theme(legend.position="none",
                            panel.background = element_rect(fill = "white", colour="black", size = 0.5, linetype = "solid")) +
                       scale_color_manual(values=c("black", "red")) + # status needs to be factor
                       geom_vline(xintercept=xintercept.neg, linetype="dotted") +
                       geom_vline(xintercept=xintercept.pos, linetype="dotted") +
                       geom_hline(yintercept=y.intercept, linetype="dotted") +
                       xlim(c(-1,1)) +
                       ylim(c(0,max.pval)) +
                       annotate("text", x=-0.75, y=3, label="DOWN", size=8) +
                       annotate("text", x=0.75, y=3, label="UP", size=8)

  if (save.as.png) {

    if (verbose) print(paste("Volcano plot: ", output.folder, "/", plot.prefix, "VolcanoPlot_Difference", Sys.Date(), ".png", sep=""))
    png(file=paste(output.folder, "/", plot.prefix, "VolcanoPlot_", plot.analysis.type, "_Difference_", Sys.Date(), ".png", sep=""), height = 600, width = 1000)
      print(volcano.plot.diff)
    dev.off()

  }

  if (display.plot) {

    #print(volcano.plot.fc)
    print(volcano.plot.diff)

  }
}
#------------------------------------------------------------------------------



