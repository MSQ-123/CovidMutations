#' Preprocess "nucmer" object to add group information
#'
#' @description Manipulate the "nucmer" object to make the analysis easier.
#'
#' @param nucmer An object called "nucmer", mutation information derived from
#' "nucmer.snp" variant file by "seqkit" software and "nucmer SNP-calling"
#' scripts.
#' @param outdir The output directory.
#' @param chinalist A list of places in China, for displacing some original cities with "China" in order
#' to make the downstream analysis easier.
#'
#' @return Saving the updated "nucmer" object.
#' @export
#' @importFrom stringr str_c
#'
#' @examples
#' data("nucmer")
#' data("chinalist")
#' #outdir <- tempdir() specify your output directory
#' nucmerr<- nucmerRMD(nucmer = nucmer, outdir = NULL, chinalist = chinalist)
nucmerRMD <- function(nucmer = nucmer, outdir = ".", chinalist = chinalist){
  nucmer <- nucmer[,c(1,2,3,4,14)]
  colnames(nucmer)<-c("rpos","rvar","qvar","qpos","ID")
  nucmer$sample <-vapply(strsplit(as.character(nucmer$ID), "[|]"), function(x) x[2], character(1))
  nucmer$time <-vapply(strsplit(as.character(nucmer$ID), "[|]"), function(x) x[3], character(1))
  nucmer$country <-vapply(strsplit(as.character(nucmer$ID), "[/]"), function(x) x[2], character(1))
  # china <-read.csv("china.txt", header = F)
  # china
  nucmer[nucmer$country %in% chinalist$V1, ]$country <-"China"
  nucmer$M_type <-str_c(nucmer$rvar,nucmer$qvar,sep ="->")
  nucmer$PM_type <-str_c(nucmer$rpos,nucmer$M_type,sep =":")
  return(nucmer)
  if(is.null(outdir) == FALSE){
    rdaPath<- file.path(outdir,"nucmerr.rda",sep = "")
    save(nucmer, file=rdaPath, compress = "xz")
  }

}
