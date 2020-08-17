#' Plot the most frequent mutational events for proteins in the SARS-CoV-2 genome
#'
#' @description Plot the most frequent mutational events for proteins selected.
#' The protein name should be specified correctly(only for SARS-CoV-2).
#'
#' @param results  The mutation effects provided by "indelSNP" function.
#' @param proteinName Proteins in the SARS-CoV-2 genome, available choices:
#' 5'UTR, NSP1~NSP10, NSP12a, NSP12b, NSP13, NSP14, NSP15, NSP16, S, ORF3a,
#' E, M, ORF6, ORF7a, ORF7b, ORF8, N, ORF10.
#' @param top The number of objects to display.
#' @param outdir The output directory.
#'
#' @return Plot the mutational events for selected proteins as output.
#' @export
#'
#' @importFrom graphics par barplot
#' @importFrom grDevices png dev.off heat.colors terrain.colors
#'
#' @examples
#' data("covid_annot")
#' covid_annot <- as.data.frame(covid_annot)
#' #outdir <- tempdir() specify your output directory
#' plotMutProteins(results = covid_annot,proteinName = "NSP2", top = 20, outdir = NULL)
plotMutProteins <- function(results = results, proteinName = "NSP2", top = 20, outdir = "."){
  covid_annot_protein <- results[results$protein == proteinName,]
  if(is.null(outdir) == TRUE){
    occ<-sort(table(covid_annot_protein$variant),decreasing=TRUE)[1:top]
    #par(las=2,mar=c(15,5,5,1))
    barplot(occ,ylab="mutation counts",main=paste("Most mutated variant", "for", proteinName,sep = " "),col=heat.colors(length(occ)))
  }
  if(is.null(outdir) == FALSE){
    figurePath<- file.path(outdir,paste("MostMutVar_","protein.png", sep = ""),sep = "")
    png(figurePath,width = 3000,height = 2000,res=300)
    occ<-sort(table(covid_annot_protein$variant),decreasing=TRUE)[1:top]
    #par(las=2,mar=c(15,5,5,1))
    barplot(occ,ylab="mutation counts",main=paste("Most mutated variant", "for", proteinName,sep = " "),col=heat.colors(length(occ)))
    dev.off()
  }

}
