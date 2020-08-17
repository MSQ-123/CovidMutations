#' Merge neighboring events of SNP, insertion and deletion.
#'
#' @description  The first step for handling the nucmer object, then effects of mutations
#' can be analysed using "indelSNP" function.
#'
#' @param nucmer An object called "nucmer", mutation information derived from
#' "nucmer.snp" variant file by "seqkit" software and "nucmer SNP-calling"
#' scripts.
#'
#' @return An updated "nucmer" object.
#' @export
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @examples
#' #The example data:
#' data("nucmer")
#' #options(stringsAsFactors = FALSE)
#'
#' #The input nucmer object can be made by the comment below:
#' #nucmer<-read.delim("nucmer.snps",as.is=TRUE,skip=4,header=FALSE)
#' #colnames(nucmer)<-c("rpos","rvar","qvar","qpos","","","","",
#' #"rlength","qlength","","","rname","qname")
#' #rownames(nucmer)<-paste0("var",1:nrow(nucmer))
#'
#' # Fix IUPAC codes
#' nucmer<-nucmer[!nucmer$qvar%in%c("B","D","H","K","M","N","R","S","V","W","Y"),]
#' nucmer<- mergeEvents(nucmer = nucmer)## This will update the nucmer object
mergeEvents <- function(nucmer = nucmer){
  samples<-unique(nucmer$qname)
  #length(samples) # 12822
  #pb<-txtProgressBar(0,length(samples),style=3)
  for (pbi in seq_along(samples)){ # This will update the nucmer object
    sample<-samples[pbi]
    allvars<-nucmer[nucmer$qname==sample,]
    snps<-allvars[(allvars[,"rvar"]!=".")&(allvars[,"qvar"]!="."),]
    inss<-allvars[(allvars[,"rvar"]=="."),]
    dels<-allvars[(allvars[,"qvar"]=="."),]
    # Merge insertions
    prevqpos<-0
    prevrowname<-NULL
    remove<-c()
    i<-1
    corrector<-0
    while(i<=nrow(inss)){
      rpos<-inss[i,"rpos"]
      rvar<-inss[i,"rvar"]
      qvar<-inss[i,"qvar"]
      qpos<-inss[i,"qpos"]
      if((qpos!=1)&(qpos==(prevqpos+1+corrector))){
        inss<-inss[-i,]
        inss[prevrowname,"qvar"]<-paste0(inss[prevrowname,"qvar"],qvar)
        corrector<-corrector+1
        i<-i-1
      } else {
        corrector<-0
        prevrowname<-rownames(inss)[i]
        prevqpos<-qpos
      }
      i<-i+1
    }
    # Merge deletions
    prevqpos<-0
    prevrowname<-NULL
    remove<-c()
    i<-1
    while(i<=nrow(dels)){
      rpos<-dels[i,"rpos"]
      rvar<-dels[i,"rvar"]
      qvar<-dels[i,"qvar"]
      qpos<-dels[i,"qpos"]

      if((qpos!=1)&(qpos==(prevqpos))){
        dels<-dels[-i,]
        dels[prevrowname,"rvar"]<-paste0(dels[prevrowname,"rvar"],rvar)
        i<-i-1
      } else {
        prevrowname<-rownames(dels)[i]
        prevqpos<-qpos
      }
      i<-i+1
    }
    # Merge SNPs
    prevqpos<-0
    prevrowname<-NULL
    remove<-c()
    i<-1
    corrector<-0
    while(i<=nrow(snps)){
      rpos<-snps[i,"rpos"]
      rvar<-snps[i,"rvar"]
      qvar<-snps[i,"qvar"]
      qpos<-snps[i,"qpos"]

      if((qpos!=1)&(qpos==(prevqpos+1+corrector))){
        snps<-snps[-i,]
        snps[prevrowname,"rvar"]<-paste0(snps[prevrowname,"rvar"],rvar)
        snps[prevrowname,"qvar"]<-paste0(snps[prevrowname,"qvar"],qvar)
        corrector<-corrector+1
        i<-i-1
      } else {
        corrector<-0
        prevrowname<-rownames(snps)[i]
        prevqpos<-qpos
      }
      i<-i+1
    }

    # Remerge back
    allvars2<-rbind(snps,inss,dels)
    remove<-setdiff(rownames(allvars),rownames(allvars2))#?setdiff
    nucmer<-nucmer[setdiff(rownames(nucmer),remove),]
    nucmer[rownames(allvars2),]<-allvars2
    #setTxtProgressBar(pb,pbi)
  }
return(nucmer)

}
