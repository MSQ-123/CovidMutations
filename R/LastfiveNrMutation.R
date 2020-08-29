#' Bacth assay analysis for last five Nr of primers
#'
#' @description Last five nucleotides of primer mutation count/type for any reverse transcription
#' polymerase chain reaction (RT-PCR) primer.
#'
#' @param nucmerr nucmerr Mutation information containing group list(derived from "nucmer" object using "nucmerRMD" function).
#' @param assays Assays dataframe including the detection ranges of mutations.
#' @param totalsample Total sample number, total cleared GISAID fasta data.
#' @param outdir The output directory. if the figurelist = TRUE, output the figure in the R session.
#' @param figurelist Whether to output the integrated plot list for each assay.
#'
#' @return Plot the mutation counts(last five nucleotides for each primer) for each assay as output.
#' @export
#' @importFrom ggplot2 aes ggplot geom_point theme_bw scale_x_continuous labs ggsave theme geom_vline element_blank
#' @importFrom cowplot plot_grid
#'
#' @examples
#' data("nucmerr")
#' data("assays")
#' totalsample <- 434
#' outdir <- tempdir()
#' LastfiveNrMutation(nucmerr = nucmerr,
#'                    assays = assays,
#'                    totalsample = totalsample,
#'                    figurelist = FALSE,
#'                    outdir = outdir)
LastfiveNrMutation <- function(nucmerr = nucmerr, assays = assays, totalsample = totalsample, figurelist = FALSE, outdir = "."){
  assays$MR_Five=0
  plist<-list()
  if(figurelist == TRUE){
    for (i in seq_along(assays$Assay)){
      #i=2
      F1=assays[i,"F1"]
      F2five=assays[i,"F2"]-5
      F2=assays[i,"F2"]
      R1=assays[i,"R1"]
      R1five=assays[i,"R1"]+5
      R2=assays[i,"R2"]

      sub_nucmer= nucmerr[nucmerr$rpos %in% c(F2five:F2,R1:R1five),]
      TMN<- length(unique(sub_nucmer$sample))
      #Total=29183 # Total Cleared GISAID fasta data
      Mutation_Ratio<- round(TMN/totalsample*100,5)

      assays[i,"MR_Five"]<- Mutation_Ratio
      rpos <- sub_nucmer$rpos
      sample<- sub_nucmer$sample
      M_type<- sub_nucmer$M_type

      plist[[i]]<- ggplot(data =sub_nucmer,aes(x=rpos,y=sample, color=M_type))+
        geom_point()+
        theme_bw()+
        scale_x_continuous(breaks=seq(F1,R2,2),limits =c(F1,R2))+
        labs(x="SARS-CoV-2 Genomic position",
             title=paste0(assays$Assay[i],"(Samples:", TMN,"/Ratio:",Mutation_Ratio,"%)"))+
        theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5),
              axis.text.y=element_blank(),
              axis.title.y=element_blank(),
              plot.title = element_text(size=8),
              legend.position="none")+
        geom_vline(aes(xintercept=F1),color="blue", linetype="dashed", size=0.5)+
        geom_vline(aes(xintercept=F2),color="blue", linetype="dashed", size=0.5)+
        geom_vline(aes(xintercept=R1),color="red", linetype="dashed", size=0.5)+
        geom_vline(aes(xintercept=R2),color="red", linetype="dashed", size=0.5)+
        geom_vline(aes(xintercept=F2five),color="gray", linetype="dashed", size=0.5)+
        geom_vline(aes(xintercept=R1five),color="gray", linetype="dashed", size=0.5)
}

    p <- plot_grid(plotlist=plist,ncol=4)
    print(p)
  }
  if(figurelist == FALSE){
    for (i in seq_along(assays$Assay)){
      #i=2
      F1=assays[i,"F1"]
      F2five=assays[i,"F2"]-5
      F2=assays[i,"F2"]
      R1=assays[i,"R1"]
      R1five=assays[i,"R1"]+5
      R2=assays[i,"R2"]

      sub_nucmer= nucmerr[nucmerr$rpos %in% c(F2five:F2,R1:R1five),]
      TMN<- length(unique(sub_nucmer$sample))
      #Total=29183 # Total Cleared GISAID fasta data
      Mutation_Ratio<- round(TMN/totalsample*100,5)

      assays[i,"MR_Five"]<- Mutation_Ratio
      rpos <- sub_nucmer$rpos
      sample<- sub_nucmer$sample
      M_type<- sub_nucmer$M_type

      p<- ggplot(data =sub_nucmer,aes(x=rpos,y=sample, color=M_type))+
        geom_point()+
        theme_bw()+
        scale_x_continuous(breaks=seq(F1,R2,2),limits =c(F1,R2))+
        labs(x="SARS-CoV-2 Genomic position",
             title=paste0(assays$Assay[i],"(Samples:", TMN,"/Ratio:",Mutation_Ratio,"%)"))+
        theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5),
              axis.text.y=element_blank(),
              axis.title.y=element_blank(),
              plot.title = element_text(size=8),
              legend.position="none")+
        geom_vline(aes(xintercept=F1),color="blue", linetype="dashed", size=0.5)+
        geom_vline(aes(xintercept=F2),color="blue", linetype="dashed", size=0.5)+
        geom_vline(aes(xintercept=R1),color="red", linetype="dashed", size=0.5)+
        geom_vline(aes(xintercept=R2),color="red", linetype="dashed", size=0.5)+
        geom_vline(aes(xintercept=F2five),color="gray", linetype="dashed", size=0.5)+
        geom_vline(aes(xintercept=R1five),color="gray", linetype="dashed", size=0.5)
      ggsave(p,filename=paste0(assays$Assay[i],'_Mutation_Count.png'),width = 12, height = 8, dpi=300,path = outdir)
    }
  }






}



