
#--------------------------------------------------------------------------------
#' Creates file name prefix used for ouptut files
#'
#' @description Returns file name prefixes as found in the annotation file (ProjectName and rRNA). File name prefix will be: RiboMethSeq_<ProjectName>_<rRNA>
#' @param annot: a data frame with general project informations. Must have columns specifying the project, rRNA, sample names, paths.
#' @param annot.project.col: a numerical value specifying the column index of the project name in the annotation <annot>. Default is 1.
#' @param annot.rna.col: a numerical value specifying the column of the rRNA in the annotation <annot>. Default is 4.
#'
#' @return file prefix used for output files (plots and .csv files)
#' @export
#'
#' @examples get.file.prefix(annotation)
get.file.prefix <- function(annot=NULL, annot.project.col=1, annot.rna.col=4){

  if (is.null(annot)) {stop("ERROR: MISSING annotation file <annot>.")}

  project.name <- unique(annot[,annot.project.col])
  if (length(project.name) > 1) stop(paste("Please check the ", annot.project.col, ".st column <Project> of your annotation file <annot>: non unique value!"))

  rRNA.name <- unique(annot[,annot.rna.col])
  if (length(rRNA.name) > 1) stop(paste("Please check the ", annot.rna.col, ".th column <rRNA> of your annotation file: non uniqe value!"))

  file.prefix <- paste("RiboMethSeq_", project.name, "_", rRNA.name, "_", sep="")

  return(file.prefix)
}
#--------------------------------------------------------------------------------




#--------------------------------------------------------------------------------
# note that the subset/select function implies the following column order which will be used for further analyses:
# (1) position, (2) type, (3) identity, (4) counts
# OTHERWISE, (if column ordering changed) you have to change the following parameters:
# <data.position.col>, <data.type.col>, <data.identity.col>, <data.counts.col>
#--------------------------------------------------------------------------------
#' Reads ribomethseq sequencing files
#'
#' @param data.folder: a character string indicating the data folder with the input files
#' @param annot: a data frame with general project informations. Must have columns specifying the project, rRNA, sample names, paths.
#' @param annot.path.col: a numerical value specifying the column with the file name in the annotation <annot>. Default is 2.
#' @param annot.name.col: a numerical value specifying the column of the sample name in the annotation <annot>. Default is 5.
#' @param annot.rna.col: a numerical value specifying the column of the rRNA in the annotation <annot>. Default is 4.
#' @param annot.group.col: a numerical value specifying the column of the biological condition in the annotation <annot>. Default is 3. Note that only 2 different biological conditions are allowed.
#' @param data.position.col: a numerical value specifying the column with the genomic position of the sequencing data in the data file. Default is 1.
#' @param data.type.col: a numerical value specifying the column with the RNA type of the sequencing data in the data file. Default is 2. If there is no column specifiying the RNA type, the parameter should set to 1.
#' @param data.identity.col: a numerical value specifying the column with the annotation of the RNA methylation site. Default is 3. If there is no column with the annotations, the parameter should set to 1.
#' @param data.counts.col: : a numerical value specifying the column of the rRNA in the annotation <annot>. Default is 4.
#' @param species: a character string specifying the species. Default is "human". In this version of the R-packge, only "human" is supported.
#' @param sep: the field separator character for the internal function read.table(). Values on each line of the file are separated by this character. If sep = "" (the default for read.table) the separator is ‘white space’, that is one or more spaces, tabs, newlines or carriage returns.
#'
#' @return List of data frames. Each element is a count table of one sample, the name of the list element corresponds to the sample name
#' @export
#'
#' @examples read.and.annotate.data(input.data.folder, project.annotation)
read.and.annotate.data <- function(data.folder=NULL, annot=NULL, species="human", annot.path.col=2, annot.group.col=3, annot.rna.col=4, annot.name.col=5, data.position.col=1, data.type.col=2, data.identity.col=3, data.counts.col=4, sep="", um14 = 15) {

  if (is.null(data.folder)) {stop("ERROR: MISSING parameter. Please provide the folder name with the input files <data.folder>.")}
  #if (is.null(input.file.format)) {stop("ERROR: MISSING parameter <input.file.format>. Should be .txt or .csv.")}
  if (is.null(annot)) {stop("ERROR: MISSING annotation file <annot>.")}
  #if (is.null(annot.path.col)) {stop("ERROR: MISSING parameter <annot.path.col>.")}
  #if (is.null(annot.name.col)) {stop("ERROR: MISSING parameter <annot.name.col>.")}

 # check annotation file
  # check RNA
  #r.rna <- unique(annot[, annot.rna.col]) # rRNA type
  if (length(unique(annot[, annot.rna.col])) > 1) {stop("ERROR: Please check your annotation file. You can only treat ONE ribosome type within an analysis!")}

  # check that there are not more than 2 biological conditions
  groups <- unique(annot[, annot.group.col]) # rRNA type
  if (length(groups) > 2) {stop("ERROR: Please check your annotation file. You can not have more than 2 biological groups!")}
  if (length(groups) == 1) {print("WARNING: Please check your annotation file. Do you really have only one biological condition?")}


  # check input folder
  lf <- lapply(annot[,annot.path.col], function(x) {
    d <- list.files(data.folder, pattern=as.character(x), recursive=T, full.names=T)
    d
  })
  lf <- unlist(lf)

  if(length(lf) == 0) {stop("Please check your input data folder <data.folder> and <annot>. The folder may not exist, may be empty, or there are no files with file names specified in the annotation.")}
  if(length(lf) != length(annot[,annot.path.col])) {stop("Please check your input data folder <data.folder> and <annot>. The annotation file does not match the file names found in the data folder.")}


  # get methylated/suspected positions
  #official.sites  <- known.sites[which(known.sites$rRNA == r.rna), c(1:3)] # select corresponding rRNA
  #suspected.sites <- suspected.sites[which(suspected.sites$rRNA == r.rna), c(1:3)]

  official.sites  <- get.known.sites.info(annot=annot, annot.rna.col=annot.rna.col, species=species) # select corresponding rRNA
  suspected.sites <- get.suspected.sites.info(annot=annot, annot.rna.col=annot.rna.col, species=species) # select corresponding rRNA

  all.samples <- lapply(lf, function(x) {
    d <- read.table(x, as.is=T, header=T, sep=sep)

    # test if new sequencing data set from CLB sequencing platform
    if (is.null(d$identity) == T){
      d <- data.frame(position=d[,data.position.col], type=NA, identity=NA, counts=d[,data.counts.col])
      # update accession variables
      data.position.col <- 1
      data.type.col <- 2
      data.identity.col <- 3
      data.counts.col <- 4
    }
    
    ### If Um14 is not in the 15 position than add to data.position.col the decalage
    if (um14 != 15) {
      n_decalage = 15 - um14
      d[,data.position.col] <- d[,data.position.col] + n_decalage
    }

    d <- subset(d, select=c(data.position.col, data.type.col, data.identity.col, data.counts.col))
    d[,2] <- unique(annot[, annot.rna.col]) # add r.rna type to data frame -> rRNA type is 2nd column (fixed)

    # IDs for known methylated sites
    if(dim(official.sites)[1] > 0){ # check that known ribosomic sites exist
      rows.official <- match(official.sites[,1], d[,1]) # select rows in df of known ribosomic sites
      d[rows.official, 3] <- stri_trim(as.character(official.sites[,3])) # replace by ID names
    }

    # IDs for suspected sites
    if(dim(suspected.sites)[1] > 0){ # check that known ribosomic sites exist
      rows.suspected <- match(suspected.sites[,1], d[,1]) # select rows in df of suspected ribosomic sites
      d[rows.suspected, 3] <- stri_trim(as.character(suspected.sites[,3])) # replace by ID names
    }
    d[,3] <- ifelse(is.na(d[,3]), d[,3], paste(d[,2], d[,3], sep="_"))
    d
  }
  )

  names(all.samples) <- get.sample.names.from.annot(lf=lf, annot=annot, annot.path.col=annot.path.col, annot.name.col=annot.name.col)

  return(all.samples)
}
#--------------------------------------------------------------------------------



#--------------------------------------------------------------------------------
#
# get sample names from annotation files
# lf =list of file paths that should match the file paths in the anntoation file / coulmn annot.path.col
# returns sample names from annotation file / column annot.name.col
#--------------------------------------------------------------------------------
#' internal function which is used by the function read.and.annotate.data()
#'
#' @param lf
#' @param annot: a data frame with general project informations. Must have columns specifying the project, rRNA, sample names, paths.
#' @param annot.path.col: a numerical value specifying the column with the file name in the annotation <annot>. Default is 2.
#' @param annot.name.col: a numerical value specifying the column of the sample name in the annotation <annot>. Default is 5.
#'
#' @return sample names
#' @export
#'
#' @examples get.sample.names.from.annot(lf, annot)
get.sample.names.from.annot <- function(lf=NULL, annot=NULL, annot.path.col=2, annot.name.col=5){

  if (is.null(annot)) {stop("ERROR: MISSING annotation file <annot>.")}
  if (is.null(lf)) {stop("ERROR: MISSING path to input folder <lf>.")}

  # check if all files are present
  check <- sapply(lf, function(x) {
    m <- tail(strsplit(x, "/")[[1]], n=1)
    m
  })
  if (! sum(check %in% annot[,annot.path.col]) == length(check)) {stop("Please check annotation file and data.folder")}

  # get sample names
  sample.names <- sapply(check, function(x){
    name <- annot[which(annot[, annot.path.col] == x), annot.name.col]
    name
  })
  return(sample.names)
}
#--------------------------------------------------------------------------------


