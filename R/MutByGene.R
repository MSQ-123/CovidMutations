#' Plot mutation counts for certain genes
#'
#' @description After annotating the mutations, this function is to plot the counts of mutational events for
#' each gene in the SARS-CoV-2 genome.
#'
#' @param nucmerr Mutation information containing group list(derived from "nucmer" object using "nucmerRMD" function).
#' @param gff3 "GFF3" format gene position data for SARS-Cov-2(the "GFF3" file should include columns named: "Gene", "Start", "Stop").
#' @param figurelist Whether to output the integrated plot list for each gene.
#' @param outdir The output directory, if the figurelist = TRUE, output the figure in the R session.
#'
#' @return Plot the mutation counts figure for each gene as output.
#' @export
#' @importFrom ggplot2 aes ggplot geom_point theme_bw scale_y_continuous labs ggsave theme
#' @importFrom cowplot plot_grid
#' @examples
#' data("nucmerr")
#' data("gene_position")
#' outdir <- tempdir()
#' MutByGene(nucmerr = nucmerr, gff3 = gene_position, figurelist = FALSE, outdir = outdir)
#' #if figurelist = TRUE, the recommendation for figure display(in pixel)is: width=1650, height=1300
MutByGene <- function(nucmerr = nucmerr, gff3 = gff3, figurelist = FALSE, outdir = "."){
  if(figurelist == FALSE){   #Requirements: gff3 includes: Gene, Start, Stop
    for ( i in seq_along(gff3$Gene)){
      P1=gff3[i,"Start"]
      P2=gff3[i,"Stop"]
      rpos <- nucmerr[nucmerr$rpos %in% P1:P2,]$rpos
      M_type <- nucmerr[nucmerr$rpos %in% P1:P2,]$M_type
      p<-ggplot(data=nucmerr[nucmerr$rpos %in% P1:P2,],aes(x=rpos, color=M_type))+
        geom_point(stat="count",size=2, alpha=2/3)+
        theme_bw()+
        scale_y_continuous(trans='log10')+
        labs(x="SARS-CoV-2 Genomic Postion",
             title = paste0(gff3$Gene[i],'_Mutation_Count'))

      ggsave(p,filename=paste0(gff3$Gene[i],'Mutation_Count.png'),width = 12, height = 8, dpi=300,path = outdir)

    }
  }
  if(figurelist == TRUE){
    ######figurelist
    # library(ggplot2)
    # library(cowplot)
    plist<-list()
    for ( i in seq_along(gff3$Gene)){
      P1=gff3[i,"Start"]
      P2=gff3[i,"Stop"]
      rpos <- nucmerr[nucmerr$rpos %in% P1:P2,]$rpos
      M_type <- nucmerr[nucmerr$rpos %in% P1:P2,]$M_type
      plist[[i]]<-ggplot(data=nucmerr[nucmerr$rpos %in% P1:P2,],aes(x=rpos, color=M_type))+
        geom_point(stat="count",size=2, alpha=2/3)+
        theme_bw()+
        scale_y_continuous(trans='log10')+
        labs(x="SARS-CoV-2 Genomic Postion",
             title = paste0(gff3$Gene[i],'_Mutation_Count'))+
        theme(legend.position="none")


    }
    plot_grid(plotlist=plist,labels = LETTERS[seq_along(gff3$Gene)])
  }
}
