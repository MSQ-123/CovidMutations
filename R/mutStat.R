#' Plot mutation statistics for nucleiotide
#'
#' @description Visualization for the top mutated samples, average mutational counts,
#' top mutated position in the genome, mutational density across the genome
#' and distribution of mutations across countries.
#'
#' @param nucmerr Mutation information containing group list(derived from "nucmer" object using "nucmerRMD" function).
#' @param outdir The output directory.
#' @param figure_Type Figure type for: "TopMuSample", "AverageMu", "TopMuPos", "MutDens", "CountryMutCount", "TopCountryMut".
#' @param type_top To plot the figure involving "top n"("TopMuSample", "TopMuPos", "TopCountryMut"), the "type_top" should specify the number of objects to display.
#' @param country To plot the figure using country as groups("CountryMutCount" and "TopCountryMut"), the "country" should be TRUE.
#' @param mutpos If the figure type is "TopCountryMut", "mutpos" can specify
#' A range of genomic position(eg. 28831:28931) for plot
#'
#' @return Plot the selected figure type as output.
#' @export
#' @importFrom ggplot2 aes ggplot geom_density theme_bw facet_grid ggsave
#' @importFrom graphics par barplot
#' @importFrom grDevices png dev.off rainbow
#' @importFrom utils head
#' @examples
#' data("nucmerr")
#' outdir <- tempdir()
#' mutStat(nucmerr = nucmerr,
#'         outdir = outdir,
#'         figure_Type = "TopCountryMut",
#'         type_top = 10,
#'         country = FALSE,
#'         mutpos = NULL)
mutStat <- function(nucmerr = nucmerr, outdir = ".", figure_Type = "TopMuSample", type_top = 10, country = FALSE, mutpos = NULL){
  # mutation statics for Nucleiotide
  if(figure_Type == "TopMuSample"){
    if(is.null(outdir) == FALSE){
      figurePath<- file.path(outdir,paste("Top", type_top, "Most_mutated_sample.png",sep = ""),sep = "")
      png(figurePath,width = 3000,height = 2000,res=300)
      ##2 Most mutated sample
      MMS<-head(sort(table(nucmerr$sample),decreasing=TRUE),n=type_top)
      par(las=3,mar=c(15,5,5,1))#par should be adjusted
      barplot(MMS,ylab="nr of mutations",main=paste("Top", type_top, "mutated samples"),col="yellow")
      dev.off()
    }
    if(is.null(outdir) == TRUE){
      MMS<-head(sort(table(nucmerr$sample),decreasing=TRUE),n=type_top)
      barplot(MMS,ylab="nr of mutations",main=paste("Top", type_top, "mutated samples"),col="yellow")
    }

  }
  ##3 Averagemutation per sample
  if(figure_Type == "AverageMu"){
    if(is.null(outdir) == FALSE){
      figurePath<- file.path(outdir,"Averagemutation_per_sample.png",sep = "")
      png(figurePath,width = 3000,height = 2000,res=300)
      MPS<- table(table(nucmerr$sample))
      par(las=1,mar=c(5,5,5,2))
      barplot(MPS, ylab="Case Number",xlab="NR of Mutation",main="Summarized Mutation Per Sample",col="lightblue")
      dev.off()
    }
    if(is.null(outdir) == TRUE){
      MPS<- table(table(nucmerr$sample))
      barplot(MPS, ylab="Case Number",xlab="NR of Mutation",main="Summarized Mutation Per Sample",col="lightblue")
    }
  }
  if(figure_Type == "TopMuPos"){
    if(is.null(outdir) == FALSE){
      figurePath<- file.path(outdir,paste("Top", type_top, "Mutation site.png", sep = ""),sep = "")
      png(figurePath,width = 3000,height = 2000,res=300)
      ##4.4 Top30 Mutation Sites
      PM_Site_top <-head(sort(table(nucmerr$PM_type), decreasing = TRUE),n=type_top)
      par(las=2,mar=c(8,5,5,2))
      barplot(PM_Site_top,ylab="number", main=paste("Top", type_top, "Mutation site"),col=rainbow(length(PM_Site_top)))
      dev.off()
    }
    if(is.null(outdir) == TRUE){
      PM_Site_top <-head(sort(table(nucmerr$PM_type), decreasing = TRUE),n=type_top)
      barplot(PM_Site_top,ylab="number", main=paste("Top", type_top, "Mutation site"),col=rainbow(length(PM_Site_top)))
    }

  }
  if(figure_Type == "MutDens"){
    ## 5 mutation density profiling by position
    MbP<-as.data.frame(table(nucmerr$rpos))
    rpos <- nucmerr$rpos
    p<- ggplot(data=nucmerr,aes(x=rpos))+
      geom_density()+
      theme_bw()
    if(is.null(outdir) == FALSE){
      ggsave(p, filename = "Mutation_Density_Bypos.png",width = 12, height = 8, dpi=300,path = outdir)
    }
    if(is.null(outdir) == TRUE){
      print(p)
    }

  }
  if(country == TRUE){
    if(figure_Type == "CountryMutCount"){
      if(is.null(outdir) == FALSE){
        figurePath<- file.path(outdir,"mutation_country_distribution.png",sep = "")
        png(figurePath,width = 3000,height = 2000,res=300)
        ### mutation country distribution

        country_sample<-nucmerr[!duplicated(nucmerr$sample),c("country", "sample")]

        #table(country_sample$country=="China")
        par(las=2,mar=c(8,5,5,2))
        barplot(sort(table(country_sample$country),decreasing =TRUE),
                col=rainbow(seq_along(country_sample$sample)))
        dev.off()
      }
      if(is.null(outdir) == TRUE){
        country_sample<-nucmerr[!duplicated(nucmerr$sample),c("country", "sample")]
        barplot(sort(table(country_sample$country),decreasing =TRUE),
                col=rainbow(seq_along(country_sample$sample)))
      }

    }
    if(figure_Type == "TopCountryMut"){
      # country mutation type in Top10 country
      if(is.null(mutpos) == TRUE){
        country_sample<-nucmerr[!duplicated(nucmerr$sample),c("country", "sample")]
        Top_country<- names(head(sort(table(country_sample$country),decreasing =TRUE),n=type_top))
        rpos<- nucmerr[nucmerr$country %in% Top_country,]$rpos
        p <- ggplot(data=nucmerr[nucmerr$country %in% Top_country,],aes(x=rpos))+
          geom_density()+
          theme_bw()+
          facet_grid(country~ .)
        if(is.null(outdir) == FALSE){
          ggsave(p, filename = paste("Top", type_top, "country_Mut.png"), width = 12, height = 8, dpi=300,path = outdir)
        }
        if(is.null(outdir) == TRUE){
          print(p)
        }

      }
      if(is.null(mutpos) == FALSE){
        country_sample<-nucmerr[!duplicated(nucmerr$sample),c("country", "sample")]
        Top_country<- names(head(sort(table(country_sample$country),decreasing =TRUE),n=type_top))
        rpos<- nucmerr[nucmerr$country %in% Top_country & nucmerr$rpos %in% mutpos,]$rpos
        p <- ggplot(data=nucmerr[nucmerr$country %in% Top_country & nucmerr$rpos %in% mutpos,],aes(x=rpos))+
          geom_density()+
          theme_bw()+
          facet_grid(country~ .)
        if(is.null(outdir) == FALSE){
          ggsave(p, filename = paste("Top", type_top, "country_MutPos.png"), width = 12, height = 8, dpi=300,path = outdir)
        }
        if(is.null(outdir) == TRUE){
          print(p)
        }

      }

    }
  }

}
