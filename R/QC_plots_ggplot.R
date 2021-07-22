#' Quality control plot
#'
#' Raw read coverage along the genome
#'
#' @param all.samples
#' @param output.folder: a character string indicating the path to output folder where the plots will be saved. Default is the working directory.
#' @param plot.prefix: a character string indicating the prefix of the file name that is used when the plot is saved as .png. Default is the empty string.
#' @param save.as.png: a logical indicating if the graphical output should be saved as .png or not. Default is TRUE.
#' @param display.plot: a logical indicating if the graphical output should be displayed at the default screen device. Default is TRUE.
#' @param data.counts.col: a numeric value specifying the index of the column of the data matrix with the counts. Default is 4.
#' @param data.position.col: a numeric value specifying the index of the column of the data matrix with the genomic position. Default is 1.
#' @param xyplot.rows
#' @param verbose: a logical indicating if the ouput messages will be printed to the console. Default is TRUE.
#'
#' @return plot
#' @export plot.quality.control.coverage
#'
#' @examples plot.quality.control.coverage(all.samples=all.samples)
plot.quality.control.coverage <- function(all.samples=NULL, output.folder="", plot.prefix="", save.as.png=TRUE, display.plot=TRUE, data.counts.col=4, data.position.col=1, xyplot.rows=2, verbose=TRUE){

  if(is.null(all.samples)) stop("ERROR: Input data <all.samples> is missing!")

  if (output.folder == "" && save.as.png) {
    output.folder <- getwd()
    if (verbose) print("No output folder provided. QC plots written in working directory.")
  }

  # define panel parameters
  xyplot.cols <- ceiling(length(all.samples)/xyplot.rows)
  #--------------------------------
  # add column with sample ID to data frame
  #--------------------------------
  df.data <- do.call("rbind", all.samples)
  sample.id <- sapply(rownames(df.data), function(x) {
    #l <- strsplit(x, "\\.")[[1]][1]
    l <- gsub("\\.[0-9]*$", "", x) # changed such that sample names may contain also points
    l
  })
  # add sample ID to data frame
  #data.sample.id.col            <- ncol(df.data) + 1
  df.data[, "sample"] <- as.factor(sample.id)

  ncol.gg <- ceiling(length(unique(sample.id)) / xyplot.rows)

  q.qc <- ggplot(data=df.data, aes(y=counts, x=position)) +
                 facet_wrap(~sample, ncol=ncol.gg) +
                 theme_bw() +
                 geom_line() +
                 xlab("Genomic position") +
                 ylab("Raw 5'/3'-end read count") +
                 theme(axis.text.x = element_text(angle=45))


  if (save.as.png) {
    if (verbose) print(paste("QC1 plot created here: ", output.folder, "/", plot.prefix, "QC_1_coverage_", Sys.Date(), ".png", sep=""))
    png(file=paste(output.folder, "/", plot.prefix, "QC_1_coverage_", Sys.Date(), ".png", sep=""))
      print(q.qc)
    dev.off()
  }
  if (display.plot) {
    print(q.qc)
  }
}




#--------------------------------
# Correlation plot
#--------------------------------
#' Quality control plot
#'
#' Plots the correlation of read counts between all samples
#'
#' @param all.samples
#' @param output.folder:  a character string indicating the path to output folder where the plots will be saved. Default is the working directory.
#' @param plot.prefix: a character string indicating the prefix of the file name that is used when the plot is saved as .png. Default is the empty string.
#' @param save.as.png: a logical indicating if the graphical output should be saved as .png or not. Default is TRUE.
#' @param display.plot: a logical indicating if the graphical output should be displayed at the default screen device. Default is TRUE.
#' @param data.counts.col: a numeric value specifying the index of the column of the data matrix with the counts. Default is 4.
#' @param data.position.col: a numeric value specifying the index of the column of the data matrix with the genomic position. Default is 1.
#' @param verbose: a logical indicating if the ouput messages will be printed to the console. Default is TRUE.
#'
#' @return plot
#' @export plot.quality.control.correlation
#'
#' @examples plot.quality.control.correlation(all.samples)
plot.quality.control.correlation <- function(all.samples=NULL, output.folder="", plot.prefix="", save.as.png=TRUE, display.plot=TRUE, data.counts.col=4, data.position.col=1, verbose=TRUE){

  if(is.null(all.samples)) stop("ERROR: Input data <df.data> is missing!")

  if (output.folder == "") {
    output.folder <- getwd()
  }

  all.samples     <- all.samples[order(names(all.samples))]
  #cor.df          <- data.frame(matrix(0, nrow=dim(all.samples[[1]])[[1]], ncol=length(all.samples) + 1))
  for (i in 1:length(all.samples)){
  all.samples[[i]] <- all.samples[[i]][,c(data.position.col,data.counts.col)]
  }
  merged.data.frame = Reduce(function(...) merge(..., all=T, by = data.position.col), all.samples)
  names(merged.data.frame)   <- c("Position", names(all.samples))
  
  ## Change the function

  #for (i in 1:length(all.samples)){
  #  out <- all.samples[[i]]
  #  if (i==1) {cor.df[,i] <- out[, data.position.col]}
  #  cor.df[,i+1 ] <- out[, data.counts.col]
  #}
  correlation <- cor(merged.data.frame[2:(dim(merged.data.frame)[1]-1),-1], method="pearson", use="complete.obs") # remove 1. and last count


  if (save.as.png) {
    if (verbose) print(paste("QC2 plot created here: ", output.folder, "/", plot.prefix, "QC_2_correlation_", Sys.Date(), ".png", sep=""))

    png(file=paste(output.folder, "/", plot.prefix, "QC_2_correlation_", Sys.Date(), ".png", sep=""))
    #opar <- par()
    corrplot(correlation, method="number", type="lower")
    #par(opar)
    dev.off()
  }
  if (display.plot) {
    #opar <- par()
    corrplot(correlation, method="number", type="lower")
    #par(opar)
  }
}
#-------------------------------



#' Quality control plots
#'
#' @param all.samples
#' @param output.folder: a character string indicating the path to output folder where the plots will be saved. Default is the working directory.
#' @param plot.prefix: a character string indicating the prefix of the file name that is used when the plot is saved as .png. Default is the empty string.
#' @param save.as.png: a logical indicating if the graphical output should be saved as .png or not. Default is TRUE.
#' @param display.plot: a logical indicating if the graphical output should be displayed at the default screen device. Default is TRUE.
#' @param data.counts.col: a numeric value specifying the index of the column of the data matrix with the counts. Default is 4.
#' @param data.position.col: a numeric value specifying the index of the column of the data matrix with the genomic position. Default is 1.
#' @param verbose: a logical indicating if the ouput messages will be printed to the console. Default is TRUE.
#'
#' @return plot
#' @export plot.QC
#'
#' @examples plot.QC(all.samples)
plot.QC <- function(all.samples=NULL, output.folder="", plot.prefix="", save.as.png=TRUE, display.plot=TRUE, data.counts.col=4, data.position.col=1, verbose=TRUE){

  plot.quality.control.coverage(all.samples       = all.samples,
                                output.folder     = output.folder,
                                plot.prefix       = plot.prefix,
                                save.as.png       = save.as.png,
                                display.plot      = display.plot,
                                data.counts.col   = data.counts.col,
                                data.position.col = data.position.col,
                                verbose           = verbose)

  plot.quality.control.correlation(all.samples       = all.samples,
                                   output.folder     = output.folder,
                                   plot.prefix       = plot.prefix,
                                   save.as.png       = save.as.png,
                                   display.plot      = display.plot,
                                   data.counts.col   = data.counts.col,
                                   data.position.col = data.position.col,
                                   verbose           = verbose)
}
