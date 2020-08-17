#' Detection of co-occurring mutations using double-assay
#'
#' @description The detection of SARS-CoV-2 is important for the
#' prevention of the outbreak and management of patients.
#' Real-time reverse-transcription polymerase chain reaction (RT-PCR)
#' assay is one of the most effective molecular diagnosis strategies
#' to detect virus in clinical laboratory. It will be more accurate
#' and practical to use double assays to detect some samples with
#' co-occurring mutations.
#'
#' @param nucmerr nucmerr Mutation information containing group list(derived from "nucmer" object using "nucmerRMD" function).
#' @param assay1 Information of the first assay(containing primers
#' locations and probe location, see the format of assays provided as
#' example data. e.g. data(assays); assay1<- assays[1,])
#' @param assay2 Information of the second assay, the format is the
#' same as the first assay.
#' @param outdir The output directory. If NULL print the plot in Rstudio.
#'
#' @return Plot three figures in a single panel, including two results of assays and a
#' "venn" plot for co-occurring mutated samples.
#' @export
#' @importFrom ggplot2 aes ggplot geom_point theme_bw scale_x_continuous geom_vline
#' theme labs scale_y_sqrt ggsave scale_y_continuous element_text geom_bar geom_text
#' @importFrom ggpubr ggarrange
#' @importFrom dplyr inner_join
#' @importFrom VennDiagram venn.diagram
#'
#' @examples
#' data("nucmerr")
#' data("assays")
#' assay1 <- assays[1,]
#' assay2 <- assays[2,]
#' #outdir <- tempdir()
#' doubleAssay(nucmerr = nucmerr,
#'             assay1 = assay1,
#'             assay2 = assay2,
#'             outdir = NULL)
doubleAssay<- function(nucmerr = nucmerr, assay1 = assay1, assay2 = assay2, outdir = "."){
  #assay1:
  F1=assay1[,"F1"]
  F2=assay1[,"F2"]
  R1=assay1[,"R1"]
  R2=assay1[,"R2"]
  P1=assay1[,"P1"]
  P2=assay1[,"P2"]

  #assay2:
  f1=assay2[,"F1"]
  f2=assay2[,"F2"]
  r1=assay2[,"R1"]
  r2=assay2[,"R2"]
  p1=assay2[,"P1"]
  p2=assay2[,"P2"]
  sub_nucmer<-subset(nucmerr, nucmerr$rpos %in% c(F1:F2,P1:P2,R1:R2))
  #write.csv(sub_nucmer,file=paste0(assays$Assay[i],'_SNP.csv'))
  sub_nucmer2 <- subset(nucmerr, nucmerr$rpos %in% c(f1:f2,p1:p2,r1:r2))

  #rm NAs
  sub_nucmer<- sub_nucmer[is.na(sub_nucmer$sample) == FALSE,]
  sub_nucmer2<- sub_nucmer2[is.na(sub_nucmer2$sample) == FALSE,]


  rpos.x <- sub_nucmer$rpos.x
  #sample<- sub_nucmer$sample
  M_type.x<- sub_nucmer$M_type.x

  overlap<- inner_join(sub_nucmer,sub_nucmer2,by = "sample")
  co_occur<- length(unique(overlap$sample))
  Mutation_Ratio1<- round(co_occur/length(unique(sub_nucmer$sample))*100,5)
  Mutation_Ratio2<- round(co_occur/length(unique(sub_nucmer2$sample))*100,5)
  assay1[,"Mutation_Ratio"]<- Mutation_Ratio1
  assay2[,"Mutation_Ratio"]<- Mutation_Ratio2
  plist <- list()
  plist[[1]]<-ggplot(data=overlap,aes(x=rpos.x, y=sample,color=M_type.x))+
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
         title=paste0(assay1$Assay,"-Co-occurring Mutant Samples:", co_occur,"/co_Mutation_Ratio:",Mutation_Ratio1,"%"))

  rpos.y <- sub_nucmer$rpos.y
  #sample<- sub_nucmer$sample
  M_type.y<- sub_nucmer$M_type.y
  plist[[2]]<-ggplot(data=overlap,aes(x=rpos.y, y=sample,color=M_type.y))+
    geom_point(size=2)+
    theme_bw()+
    scale_x_continuous(breaks=seq(f1,r2,2),limits =c(f1,r2))+
    geom_vline(aes(xintercept=F1),color="blue", linetype="dashed", size=0.5)+
    geom_vline(aes(xintercept=F2),color="blue", linetype="dashed", size=0.5)+
    geom_vline(aes(xintercept=R1),color="red", linetype="dashed", size=0.5)+
    geom_vline(aes(xintercept=R2),color="red", linetype="dashed", size=0.5)+
    geom_vline(aes(xintercept=P1),color="gray", linetype="solid", size=0.5)+
    geom_vline(aes(xintercept=P2),color="gray", linetype="solid", size=0.5)+
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))+
    labs(x="SARS-CoV-2 Genomic position",
         title=paste0(assay2$Assay,"-Co-occurring Mutant Samples:", co_occur,"/co_Mutation_Ratio:",Mutation_Ratio2,"%"))



  plist[[3]] <- venn.diagram(
    list(assay1 = unique(sub_nucmer$sample),  #区域1的数
         assay2 = unique(sub_nucmer2$sample)),   #区域2的数
    #cross.area = unique(overlap$sample)),
    resolution = 300, imagetype = "png", alpha=c(0.5,0.5),
    fill=c("red","yellow"), cat.fontface=4,
    main="Intersection of identified mutations by different assays",
    main.cex = 2, #main.fontfamily = "Times New Roman",
    filename = NULL
  )
  p<- ggarrange(plist[[1]],plist[[2]],plist[[3]],ncol=2,nrow=2,labels=c("A","B","C"))
  if(is.null(outdir) == TRUE){
    print(p)
  }
  if(is.null(outdir) == FALSE){
    ggsave(p,filename=paste0("doubleAssay",'.png'),width = 18, height = 12, dpi=300,path = outdir)
    write.csv(overlap, file.path(outdir, "overlap_samples.csv"))
  }
}
