#' Global mutational events profiling of proteins
#'
#' @description This function is to visualize the global protein mutational pattern in the
#' SARS-CoV-2 genome.
#'
#' @param covid_annot The mutation effects provided by "indelSNP" function.
#' @param outdir The output directory.
#' @param figure_Type Figure type for either "heatmap" or "count".
#' @param top The number of variants to plot.
#' @param country Choose a country to plot the mutational pattern or
#' choose "global" to profile mutations across all countries. The
#' default is "global".
#'
#' @return Plot the selected figure type as output.
#' @export
#' @importFrom ggplot2 aes ggplot geom_point theme_bw labs scale_y_sqrt ggsave scale_y_continuous
#'
#' @examples
#' data("covid_annot")
#' outdir <- tempdir()
#' # make sure the covid_annot is a dataframe
#' covid_annot <- as.data.frame(covid_annot)
#' globalProteinMut(covid_annot = covid_annot,
#'                  outdir = outdir,
#'                  figure_Type = "heatmap",
#'                  top = 10,
#'                  country = "USA")
globalProteinMut <- function(covid_annot = covid_annot, outdir = ".", figure_Type = "heatmap", top = 10, country = "global"){
  covid_annot$country <-vapply(strsplit(as.character(covid_annot$sample), "[/]"), function(x) x[2], character(1))
  if(country != "global"){
    covid_annot<-covid_annot[covid_annot$country == country, ]
    covid_annot<- covid_annot[covid_annot$variant %in% names(head(sort(table(covid_annot$variant),decreasing=TRUE),n=top)),]
    if(figure_Type == "heatmap"){
      qpos <- covid_annot$qpos
      sample <- covid_annot$sample
      variant <- covid_annot$variant
      p <- ggplot(data=covid_annot,aes(x=qpos, y=vapply(strsplit(as.character(sample), "[|]"), function(x) x[2], character(1))))+
        geom_point(size=0.001, alpha=2/3,aes(color=variant))+
        theme_bw()+
        labs(x="SARS-CoV-2 Genomic position")
      if(is.null(outdir) == FALSE){
        ggsave(p,filename = paste("Top", top, country, "_protein.png", sep = ""), width = 12, height = 8, dpi=300,path = outdir)
      }
      if(is.null(outdir) == TRUE){
        print(p)
      }

    }
    if(figure_Type == "count"){
      qpos <- covid_annot$qpos
      sample <- covid_annot$sample
      variant <- covid_annot$variant
      p<-ggplot(data=covid_annot,aes(x=qpos,color=variant))+
        geom_point(stat="count",size=2, alpha=2/3)+
        theme_bw()+
        labs(x="SARS-CoV-2 Genomic position")+
        scale_y_continuous(trans='log10')+
        labs(x="SARS-CoV-2 Genomic Postion",
             title = 'SARS-CoV-2: Mutation_Count')
      if(is.null(outdir) == FALSE){
        ggsave(p,filename = paste("Top", 10, country, "_protein_count.png", sep = ""),width = 12, height = 8, dpi=300, path = outdir)
      }
      if(is.null(outdir) == TRUE){
        print(p)
      }

    }
  }
  if(country == "global"){
    covid_annot<- covid_annot[covid_annot$variant %in% names(head(sort(table(covid_annot$variant),decreasing=TRUE),n=top)),]
    if(figure_Type == "heatmap"){
      qpos <- covid_annot$qpos
      sample <- covid_annot$sample
      variant <- covid_annot$variant
      p <- ggplot(data=covid_annot,aes(x=qpos, y=vapply(strsplit(as.character(sample), "[|]"), function(x) x[2], character(1))))+
        geom_point(size=0.001, alpha=2/3,aes(color=variant))+
        theme_bw()+
        labs(x="SARS-CoV-2 Genomic position")
      if(is.null(outdir) == FALSE){
        ggsave(p,filename = paste("Top", top, "_Global_protein.png", sep = ""), width = 12, height = 8, dpi=300,path = outdir)
      }
      if(is.null(outdir) == TRUE){
        print(p)
      }

    }
    if(figure_Type == "count"){
      qpos <- covid_annot$qpos
      sample <- covid_annot$sample
      variant <- covid_annot$variant
      p<-ggplot(data=covid_annot,aes(x=qpos,color=variant))+
        geom_point(stat="count",size=2, alpha=2/3)+
        theme_bw()+
        labs(x="SARS-CoV-2 Genomic position")+
        scale_y_continuous(trans='log10')+
        labs(x="SARS-CoV-2 Genomic Postion",
             title = 'SARS-CoV-2: Mutation_Count')
      if(is.null(outdir) == FALSE){
        ggsave(p,filename = paste("Top", 10, "_Global_protein_count.png", sep = ""),width = 12, height = 8, dpi=300, path = outdir)
      }
      if(is.null(outdir) == TRUE){
        print(p)
      }

    }
  }

}
