#' Calculate the mutation detection rate using different assays
#'
#' @description This function is to use the well established assays information to detect mutations in different
#' SARS-CoV-2 genomic sites. The output will be series of figures presenting the mutation profile using a
#' specific assay and a figure for comparison between the mutation detection rate in each primers binding region.
#'
#' @param nucmerr nucmerr Mutation information containing group list(derived from "nucmer" object using "nucmerRMD" function).
#' @param assays Assays dataframe including the detection ranges of mutations.
#' @param totalsample Total sample number, total cleared GISAID fasta data.
#' @param plotType Figure type for either "barplot" or "logtrans".
#' @param outdir The output directory.
#'
#' @return Plot the selected figure type as output.
#' @export
#' @importFrom ggplot2 aes ggplot geom_point theme_bw scale_x_continuous geom_vline
#' theme labs scale_y_sqrt ggsave scale_y_continuous element_text geom_bar geom_text
#' annotate arrow
#' @importFrom graphics par barplot
#' @importFrom grDevices png dev.off heat.colors rainbow
#'
#'
#' @examples
#' data("nucmerr")
#' data("assays")
#' Total <- 434 ## Total Cleared GISAID fasta data, sekitseq
#' outdir <- tempdir()
#' #Output the results
#' AssayMutRatio(nucmerr = nucmerr,
#'               assays = assays,
#'               totalsample = Total,
#'               plotType = "logtrans",
#'               outdir = outdir)

AssayMutRatio <- function(nucmerr = nucmerr, assays = assays, totalsample = totalsample, plotType = "barplot", outdir = "."){
  for (i in seq_along(assays$Assay)){
    F1=assays[i,"F1"]
    F2=assays[i,"F2"]
    R1=assays[i,"R1"]
    R2=assays[i,"R2"]
    P1=assays[i,"P1"]
    P2=assays[i,"P2"]

    sub_nucmer<-subset(nucmerr, nucmerr$rpos %in% c(F1:F2,P1:P2,R1:R2))
    #write.csv(sub_nucmer,file=paste0(assays$Assay[i],'_SNP.csv'))

    TMN<- length(unique(sub_nucmer$sample))
    #Total=29183  # Total Cleared GISAID fasta data, sekitseq
    #grep -c '^>' Gisaid_RMD.fasta  Total = total fasta seq
    #Mutation_Ratio = assay_identified_MutSample/total sample number
    #Total=37527
    Mutation_Ratio<- round(TMN/totalsample*100,5)
    assays[i,"Mutation_Ratio"]<- Mutation_Ratio

    rpos <- sub_nucmer$rpos
    sample<- sub_nucmer$sample
    M_type<- sub_nucmer$M_type

    p<-ggplot(data=sub_nucmer,aes(x=rpos, y=sample,color=M_type))+
      geom_point(size=2)+
      theme_bw()+
      scale_x_continuous(breaks=seq(F1,R2,2),limits =c(F1,R2))+
      geom_vline(aes(xintercept=F1),color="blue", linetype="dashed", size=0.5)+
      geom_vline(aes(xintercept=F2),color="blue", linetype="dashed", size=0.5)+
      geom_vline(aes(xintercept=R1),color="red", linetype="dashed", size=0.5)+
      geom_vline(aes(xintercept=R2),color="red", linetype="dashed", size=0.5)+
      geom_vline(aes(xintercept=P1),color="gray", linetype="solid", size=0.5)+
      geom_vline(aes(xintercept=P2),color="gray", linetype="solid", size=0.5)+
      theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))+
      labs(x="SARS-CoV-2 Genomic position",
           title=paste0(assays$Assay[i],"-Total Mutant Samples:", TMN,"/Mutation_Ratio:",Mutation_Ratio,"%"))
    P<- p+annotate("segment",x=F1,xend=F2,y=0,yend = 0, arrow = arrow(ends = "last", type = "closed"), colour = "red", size = 1)+annotate("segment",x=R2,xend=R1,y=0,yend = 0, arrow = arrow(ends = "last", type = "closed"), colour = "red", size = 1)+annotate("segment",x=P1,xend=P2,y=0,yend = 0, arrow = arrow(ends = "last", type = "closed"), colour="blue", size = 1)
    if(is.null(outdir) == FALSE){
      ggsave(P,filename=paste0(assays$Assay[i],'.png'),width = 12, height = 8, dpi=300,path = outdir)
    }
  }




  ######
  # csvPath<- file.path(outdir,"Assays#2.csv",sep = "")
  # write.csv(assays, file=csvPath, row.names = F,)
  if(plotType == "barplot"){
    if(is.null(outdir) == FALSE){
      figurePath<- file.path(outdir,"Mutation_Ratio_barplot.png",sep = "")
      png(figurePath,width = 3000,height = 2000,res=300)
      #par(las=3,mar=c(8,5,5,2))
      barplot(data=assays, log2(assays$Mutation_Ratio) ~ assays$Assay, col=heat.colors(12))
      dev.off()
    }
    if(is.null(outdir) == TRUE){
      barplot(data=assays, log2(assays$Mutation_Ratio) ~ assays$Assay, col=heat.colors(12))
    }

  }
  if(plotType == "logtrans"){
    # library(ggsci)

    order_assay=assays[order(assays$Mutation_Ratio),]$Assay

    assays$Assay=factor(assays$Assay, levels=order_assay)
    p <- ggplot(data=assays,aes(x=assays$Assay,y=assays$Mutation_Ratio,fill=assays$Assay))+
      geom_bar(stat="identity")+
      geom_text(aes(label=assays$Mutation_Ratio),vjust=0,angle = 90)+
      theme_bw()+
      theme(axis.text.x = element_text(angle = 45,hjust = 1,size=12, face="bold"),
            axis.text.y = element_text(size=12,face="bold"))+
      scale_y_sqrt()
    scale_y_continuous(trans='log10')
    if(is.null(outdir) == FALSE){
      ggsave(p,filename="log2_Mutation_Ratio.png",width = 12, height = 8, dpi=300,path = outdir)
    }
    if(is.null(outdir) == TRUE){
      print(p)
    }

  }

}
