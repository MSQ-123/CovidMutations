#' Provide effects of each single nucleotide polymorphism (SNP), insertion and deletion in virus genome
#'
#' @description This function is to annotate the mutational events and indicate their potential effects on the proteins.
#' Mutational events include SNP, insertion and deletion.
#'
#' @param nucmer An object called "nucmer", mutation information derived
#' from "nucmer.snp" variant file by "seqkit" software and "nucmer
#' SNP-calling" scripts. To be processed by "indelSNP" function, The nucmer
#' object should be first transformed by "mergeEvents" function.
#' @param saveRda Whether to save the results as ".rda" file.
#' @param refseq SARS-Cov-2 genomic reference sequence.
#' @param gff3 "GFF3" format annotation data for SARS-Cov-2.
#' @param annot Annotation of genes(corresponding proteins) list from "GFF3" file by "setNames(gff3[,10],gff3[,9])".
#' @param outdir The output directory.
#'
#' @return Write the result as ".csv" file to the specified directory.
#' @export
#' @importFrom seqinr read.fasta translate
#' @importFrom utils txtProgressBar setTxtProgressBar write.csv
#' @examples
#' data("nucmer")
#' # Fix IUPAC codes
#' nucmer<-nucmer[!nucmer$qvar%in%c("B","D","H","K","M","N","R","S","V","W","Y"),]
#' nucmer<- mergeEvents(nucmer = nucmer)## This will update the nucmer object
#' data("refseq")
#' data("gff3")
#' annot <- setNames(gff3[,10],gff3[,9])
#' outdir <- tempdir()
#' indelSNP(nucmer = nucmer,
#'          saveRda = FALSE,
#'          refseq = refseq,
#'          gff3 = gff3,
#'          annot = annot,
#'          outdir = outdir)
indelSNP <- function(nucmer = nucmer, saveRda = FALSE, refseq = refseq, gff3 = gff3, annot = annot, outdir = "."){
  header<-c("sample","refpos","refvar","qvar","qpos","qlength","protein","variant","varclass","annotation")
  results<-matrix(NA,ncol=length(header),nrow=0)
  colnames(results)<-header
  samples<-unique(nucmer$qname)
  # cod <- c("GCT","GCC","GCA","GCG","TGT","TGC","GAT","GAC","GAA","GAG","TTT",
  #          "TTC","GGT","GGC","GGA","GGG","CAT","CAC","ATA","ATT","ATC","AAA",
  #          "AAG","TTA","TTG","CTT","CTC","CTA","CTG","ATG","AAT","AAC","CCT",
  #          "CCC","CCA","CCG","CAA","CAG","CGT","CGC","CGA","CGG","AGA","AGG",
  #          "TCT","TCC","TCA","TCG","AGT","AGC","ACT","ACC","ACA","ACG","GTT",
  #          "GTC","GTA","GTG","TGG","TAT","TAC","TAA","TAG","TGA")
  # aa <- c("A","A","A","A","C","C","D","D","E","E","F","F","G","G","G","G","H",
  #         "H","I","I","I","K","K","L","L","L","L","L","L","M","N","N","P","P",
  #         "P","P","Q","Q","R","R","R","R","R","R","S","S","S","S","S","S","T",
  #         "T","T","T","V","V","V","V","W","Y","Y","*","*","*")
  # codon <- data.frame(cod,aa)
  #pb<-txtProgressBar(0,length(samples),style=3)
  for (pbi in 1:length(samples)){ # This will update the nucmer object
    #add condon list, the stop codon is for "*".

    sample<-samples[pbi]
    allvars<-nucmer[nucmer$qname==sample,]
    # Check changes in query protein sequence according to variants
    for(i in 1:nrow(allvars)){ # Assuming they are sorted numerically
      nucline<-allvars[i,]
      rpos<-nucline[1,"rpos"]
      rvar<-nucline[1,"rvar"]
      qvar<-nucline[1,"qvar"]
      qpos<-nucline[1,"qpos"]
      qlength<-nucline[1,"qlength"]

      # Match over GFF3 annotation   for each var in a sample: match the var with the GFF3
      a<-rpos-gff3[,4]
      b<-rpos-gff3[,5]
      signs<-sign(a)*sign(b)
      w<-which(signs==-1)#if signs == -1 then  rpos in GFF anno

      # Outside genes scenarios
      if(length(w)==0){
        if(rpos<gff3[1,4]){    #rpos in the upstream of all gene
          protein<-"5'UTR";output<-c(rpos,"extragenic")
        } else if(rpos>gff3[nrow(gff3),5]){  #?rpos in the downstream of all gene rpos>gff3[nrow(gff3),5]
          protein<-"3'UTR";output<-c(rpos,"extragenic")
        } else {
          protein<-"intergenic";output<-c(rpos,"extragenic")
        }

      } else{ # Inside genes scenario
        start<-gff3[w,4]
        end<-gff3[w,5]
        protein<-gff3[w,9]


        #refdnaseq<-DNAString(paste0(refseq[start:end],collapse=""))
        refdnaseq <- paste0(refseq[start:end],collapse="")
        #seqm<- toupper(refdnaseq)

        refpepseq<- translate(strsplit(refdnaseq, split= "")[[1]])
        #refpepseq<-Biostrings::translate(refdnaseq) #obtain protein seq
        #refpepseq<-strsplit(as.character(refpepseq),"")[[1]]
        if(qvar=="."){ # Deletion scenario
          if((nchar(rvar)%%3)!=0){ # Deletion frameshift scenario
            mutpos<-round((rpos-start+1)/3)
            output<-c(paste0(refpepseq[mutpos],mutpos),"deletion_frameshift")
          } else { # In-frame deletion
            varseq<-refseq
            varseq<-varseq[-(rpos:(rpos+nchar(rvar)-1))]
            varseq<-varseq[start:(end-nchar(rvar))]
            #vardnaseq<-DNAString(paste0(varseq,collapse=""))
            vardnaseq <- paste0(varseq,collapse="")

            varpepseq<- translate(strsplit(vardnaseq, split= "")[[1]])
            #varpepseq<-Biostrings::translate(vardnaseq)
            #varpepseq<-strsplit(as.character(varpepseq),"")[[1]]

            for(j in 1:length(refpepseq)){
              refj<-refpepseq[j]#
              varj<-varpepseq[j]
              if(refj!=varj){
                if(varj=="*"){
                  output<-c(paste0(refj,j),"deletion_stop")
                } else {
                  output<-c(paste0(refj,j),"deletion")
                }
                break()
              }
            }
          }
        } else if(rvar=="."){ # Insertion scenario
          if((nchar(qvar)%%3)!=0){ # Insertion frameshift scenario
            mutpos<-round((rpos-start+1)/3)
            output<-c(paste0(refpepseq[mutpos],mutpos),"insertion_frameshift")
          } else { # In-frame insertion
            varseq<-c(refseq[1:rpos],strsplit(qvar,"")[[1]],refseq[(rpos+1):length(refseq)])
            varseq<-varseq[start:(end+nchar(qvar))]
            #vardnaseq<-DNAString(paste0(varseq,collapse=""))

            vardnaseq <- paste0(varseq,collapse="")

            varpepseq<- translate(strsplit(vardnaseq, split= "")[[1]])
            #varpepseq<-Biostrings::translate(vardnaseq)
            #varpepseq<-strsplit(as.character(varpepseq),"")[[1]]

            for(j in 1:length(refpepseq)){
              refj<-refpepseq[j]
              varj<-varpepseq[j]
              if(refj!=varj){
                nr_aa_inserted<-nchar(qvar)/3
                multivarj<-varpepseq[j:(j+nr_aa_inserted-1)]
                if(any(multivarj=="*")){
                  multivarj<-paste0(multivarj,collapse="")
                  output<-c(paste0(multivarj,j),"insertion_stop")
                } else{
                  multivarj<-paste0(multivarj,collapse="")
                  output<-c(paste0(multivarj,j),"insertion")
                }
                break()
              }
            }
          }
        } else { # SNP scenario
          if(nchar(qvar)==1){ # ?Single nucleotide scenario
            varseq<-refseq
            varseq[rpos]<-qvar
            varseq<-varseq[start:end]
            #vardnaseq<-DNAString(paste0(varseq,collapse=""))

            vardnaseq <- paste0(varseq,collapse="")

            varpepseq<- translate(strsplit(vardnaseq, split= "")[[1]])
            #varpepseq<-Biostrings::translate(vardnaseq)
            #varpepseq<-strsplit(as.character(varpepseq),"")[[1]]
            mutpos<-which(varpepseq!=refpepseq)  #split
            if(length(mutpos)==0){ # Silent SNP scenario
              mutpos<-round((rpos-start+1)/3)
              refaa<-refpepseq[mutpos]
              varaa<-varpepseq[mutpos]
              output<-c(paste0(refaa,mutpos,varaa),"SNP_silent")
            } else { # Changed aa scenario
              refaa<-refpepseq[mutpos]
              varaa<-varpepseq[mutpos]
              if(varaa=="*"){
                output<-c(paste0(refaa,mutpos,varaa),"SNP_stop")
              } else {
                output<-c(paste0(refaa,mutpos,varaa),"SNP")
              }
            }
          } else { # Multiple neighboring nucleotides
            varlength<-nchar(qvar)
            varseq<-refseq
            varseq[rpos:(rpos+varlength-1)]<-strsplit(qvar,"")[[1]]
            varseq<-varseq[start:end]
            #vardnaseq<-DNAString(paste0(varseq,collapse=""))

            vardnaseq <- paste0(varseq,collapse="")

            varpepseq<- translate(strsplit(vardnaseq, split= "")[[1]])
            #varpepseq<-Biostrings::translate(vardnaseq)
            #varpepseq<-strsplit(as.character(varpepseq),"")[[1]]
            mutpos<-which(varpepseq!=refpepseq)
            if(length(mutpos)==0){ # Silent SNP scenario
              mutpos<-round((rpos-start+1)/3)
              refaa<-refpepseq[mutpos]
              varaa<-varpepseq[mutpos]
              output<-c(paste0(refaa,mutpos,varaa),"SNP_silent")
            } else { # Changed aa scenario
              refaa<-paste0(refpepseq[mutpos],collapse="")
              varaa<-paste0(varpepseq[mutpos],collapse="")
              if(any(varaa=="*")){
                output<-c(paste0(refaa,mutpos[1],varaa),"SNP_stop")
              } else {
                output<-c(paste0(refaa,mutpos[1],varaa),"SNP")
              }
            }
          }
        }
      }
      results<-rbind(results,c(sample,rpos,rvar,qvar,qpos,qlength,protein,output,annot[protein]))
    }
    #setTxtProgressBar(pb,pbi)
  }
  results<-as.data.frame(results,stringsAsFactors=FALSE)
  return(results)
  csvPath<- file.path(outdir,"covid_annot.csv",sep = "")
  write.csv(results,file=csvPath,row.names = FALSE)
  if(saveRda == TRUE){
    rdaPath <- file.path(outdir,"covid_annot.rda",sep = "")
    save(results,file=rdaPath)
  }
}
