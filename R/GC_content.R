
#-------------------------------
# load GC content for human ribosome
#-------------------------------
#' Load GC content and positions of methylated sites for human ribosomes
#'
#' @param annot
#' @param annot.rna.col
#'
#' @return List of GC content of known methylated and suspected positions
#' @export
#'
#' @examples load.gc.content.for.human.ribosome(annot)
load.gc.content.for.human.ribosome <- function(annot=NULL, annot.rna.col=4){

  if (is.null(annot)) {stop("MISSING parameter: please provide project annotation file <annot>!")}

  ribosome <- unique(annot[, annot.rna.col])
  if (length(ribosome) > 1) {stop("Please verify your annotation file <annot> and the rRNA column. You can only treat one ribosome type within an analysis!")}

  if (! toupper(ribosome) %in% c("28S", "5.8S", "18S", "5S")) {stop("Ribosome type ", ribo, " has to be one of: 5S, 5.8S, 18S, 28S")}

  if (toupper(ribosome) == "28S"){

      gc.probas       <- wdw.probas.28s
      pos.meth.vector <- human.methylated[which(human.methylated$rRNA == "28S"), 1]
      pos.suspected   <- human.suspected[which(human.suspected$rRNA == "28S"), 1]

  } else if(toupper(ribosome) == "5S"){

      gc.probas       <- wdw.probas.5s
      pos.meth.vector <- c(0)
      pos.suspected   <- c(0)
      print("No known methylated and no suspected sites for ribosomic protein 5S.")

  } else if(toupper(ribosome) == "5.8S"){

      gc.probas       <- wdw.probas.5.8s
      pos.meth.vector <- human.methylated[which(human.methylated$rRNA == "5.8S"), 1]
      pos.suspected   <- c(0)
      print("No suspected sites for ribosomic protein 5.8S.")

  } else if(toupper(ribosome) == "18S"){

      gc.probas       <- wdw.probas.18s
      pos.meth.vector <- human.methylated[which(human.methylated$rRNA == "18S"), 1]
      pos.suspected   <- human.suspected[which(human.suspected$rRNA == "18S"), 1]
  }

  return(list(gc.probas=gc.probas, pos.meth.vector=pos.meth.vector, pos.suspected=pos.suspected))
}
#----------------------


#' Load GC content and positions of methylated sites for human ribosomes
#'
#' @param annot
#' @param annot.rna.col
#'
#' @return List of GC content of known methylated and suspected positions
#' @export
#'
#' @examples load.gc.content.for.mouse.ribosome(annot)
load.gc.content.for.mouse.ribosome <- function(annot=NULL, annot.rna.col=4){

  if (is.null(annot)) {stop("MISSING parameter: please provide project annotation file <annot>!")}

  ribosome <- unique(annot[, annot.rna.col])
  if (length(ribosome) > 1) {stop("Please verify your annotation file <annot> and the rRNA column. You can only treat one ribosome type within an analysis!")}

  if (! toupper(ribosome) %in% c("28S", "5.8S", "18S")) {stop("Ribosome type ", ribo, " has to be one of: 5.8S, 18S, 28S")}

  if (toupper(ribosome) == "28S"){

    gc.probas       <- mouse.wdw.probas.28s
    pos.meth.vector <- mouse.methylated[which(mouse.methylated$rRNA == "28S"), 1]
    pos.suspected   <- mouse.suspected[which(mouse.suspected$rRNA == "28S"), 1]

  } else if(toupper(ribosome) == "5.8S"){

    gc.probas       <- mouse.wdw.probas.5.8s
    pos.meth.vector <- mouse.methylated[which(mouse.methylated$rRNA == "5.8S"), 1]
    pos.suspected   <- mouse.suspected[which(mouse.suspected$rRNA == "5.8S"), 1]
    print("No suspected sites for ribosomic protein 5.8S.")

  } else if(toupper(ribosome) == "18S"){

    gc.probas       <- mouse.wdw.probas.18s
    pos.meth.vector <- mouse.methylated[which(mouse.methylated$rRNA == "18S"), 1]
    pos.suspected   <- mouse.suspected[which(mouse.suspected$rRNA == "18S"), 1]
  }

  return(list(gc.probas=gc.probas, pos.meth.vector=pos.meth.vector, pos.suspected=pos.suspected))
}
#----------------------






#----------------------
#' Get position vector of methylated positions (internal function)
#'
#' @param annot
#' @param annot.rna.col
#' @param species
#'
#' @return Returns a vector of known methylated positions (if existing) of the ribosomic unit
#' @export
#'
#' @examples get.known.sites.vec(annot)
get.known.sites.info <- function(annot=NULL, annot.rna.col=4, species="human"){

  if (is.null(annot)) {stop("MISSING parameter (within get.known.sites.vec()): please provide project annotation file <annot>!")}

  if (tolower(species) == "human"){

    sites <- human.methylated[which(human.methylated[, "rRNA"] == as.character(unique(annot[,annot.rna.col]))), c(1:3)]  # 2nd list element, corresponding to pos.meth.vector

  } else if (tolower(species) == "mouse"){

    sites <- mouse.methylated[which(mouse.methylated[, "rRNA"] == as.character(unique(annot[,annot.rna.col]))), c(1:3)]  # 2nd list element, corresponding to pos.meth.vector

  } else {

      print("ERROR: Please use a valid species. Choose among: <human> or <mouse>. Other species under development :-)")
  }

  return(sites)
}
#----------------------

#----------------------
#' Get vector of positions suspected to be methylated (internal function)
#'
#' @param annot
#' @param annot.rna.col
#' @param species
#'
#' @return Returns positions suspected to be methylated (if existing) of the ribosomic unit
#' @export
#'
#' @examples get.suspected.sites.info(annot)
get.suspected.sites.info <- function(annot=NULL, annot.rna.col=4, species="human"){

  if (is.null(annot)) {stop("MISSING parameter (within get.suspected.sites.vec()): please provide project annotation file <annot>!")}

  if (tolower(species) == "human"){

    sites <- human.suspected[which(human.suspected[, "rRNA"] == as.character(unique(annot[,annot.rna.col]))), c(1:3)] # 3rd list element, corresponding to pos.suspected

  } else if (tolower(species) == "mouse"){

    sites <- mouse.suspected[which(mouse.suspected[, "rRNA"] == as.character(unique(annot[,annot.rna.col]))), c(1:3)] # 3rd list element, corresponding to pos.suspected

  } else {

    print("ERROR: Please use a valid species. Choose among: <human>. Other species under development :-)")
  }

  return(sites)
}
#----------------------
