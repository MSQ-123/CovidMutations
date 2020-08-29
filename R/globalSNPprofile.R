#' Global single nucleotide polymorphism (SNP) profiling in virus genome
#'
#' @description This function is to visualize the global SNP pattern in the
#' SARS-CoV-2 genome.
#'
#' @param nucmerr Mutation information containing group list(derived from "nucmer" object using "nucmerRMD" function).
#' @param outdir The output directory.
#' @param figure_Type Figure type for either "heatmap" or "count".
#' @param country Choose a country to plot the mutational pattern or
#' choose "global" to profile mutations across all countries. The
#' default is "global".
#' @param top The number of mutational classes to plot.
#'
#' @return Plot the selected figure type as output.
#' @export
#' @importFrom ggplot2 aes ggplot geom_point theme_bw labs scale_y_sqrt ggsave scale_y_continuous
#'
#' @examples
#' data("nucmerr")
#' outdir <- tempdir()
#' globalSNPprofile(nucmerr = nucmerr,
#'                  outdir = outdir,
#'                  figure_Type = "heatmap",
#'                  country = "global",
#'                  top = 5)
globalSNPprofile <- function(nucmerr = nucmerr, outdir = ".", figure_Type = "heatmap", country = "global", top = 5){
  if(country != "global"){
    sub_nucmer<-nucmerr[nucmerr$country == country, ]
    rpos <- sub_nucmer$rpos
    sample <- sub_nucmer$sample
    M_type <- sub_nucmer$M_type
    sub_nucmer<- sub_nucmer[sub_nucmer$M_type %in% names(head(sort(table(sub_nucmer$M_type),decreasing=TRUE),n=top)),]
    if(figure_Type == "heatmap"){

      p<-ggplot(data=sub_nucmer,aes(x=rpos, y=sample))+
        geom_point(size=0.001, alpha=2/3,aes(color=M_type))+
        theme_bw()+
        labs(x="SARS-CoV-2 Genomic position")
      if(is.null(outdir) == FALSE){
        ggsave(p,filename = paste(country, "_SNP.png", sep = ""), width = 12, height = 8, dpi=300,path = outdir)
      }
      if(is.null(outdir) == TRUE){
        print(p)
      }

    }
    if(figure_Type == "count"){
      p<-ggplot(data=sub_nucmer,aes(x=rpos,color=M_type))+
        geom_point(stat="count",size=2, alpha=2/3)+
        theme_bw()+
        labs(x="SARS-CoV-2 Genomic position")+
        scale_y_continuous(trans='log10')+
        labs(x="SARS-CoV-2 Genomic Postion",
             title = 'SARS-CoV-2: Mutation_Count')
      if(is.null(outdir) == FALSE){
        ggsave(p,filename = paste(country, "_SNP_count.png", sep = ""),width = 12, height = 8, dpi=300, path = outdir)
      }
      if(is.null(outdir) == TRUE){
        print(p)
      }

    }
  }
  if(country == "global"){
    sub_nucmer <- nucmerr
    rpos <- sub_nucmer$rpos
    sample <- sub_nucmer$sample
    M_type <- sub_nucmer$M_type
    sub_nucmer<- sub_nucmer[sub_nucmer$M_type %in% names(head(sort(table(sub_nucmer$M_type),decreasing=TRUE),n=top)),]
    if(figure_Type == "heatmap"){

      p<-ggplot(data=sub_nucmer,aes(x=rpos, y=sample))+
        geom_point(size=0.001, alpha=2/3,aes(color=M_type))+
        theme_bw()+
        labs(x="SARS-CoV-2 Genomic position")
      if(is.null(outdir) == FALSE){
        ggsave(p,filename = paste(country, "_SNP.png", sep = ""), width = 12, height = 8, dpi=300,path = outdir)
      }
      if(is.null(outdir) == TRUE){
        print(p)
      }

    }
    if(figure_Type == "count"){
      p<-ggplot(data=sub_nucmer,aes(x=rpos,color=M_type))+
        geom_point(stat="count",size=2, alpha=2/3)+
        theme_bw()+
        labs(x="SARS-CoV-2 Genomic position")+
        scale_y_continuous(trans='log10')+
        labs(x="SARS-CoV-2 Genomic Postion",
             title = 'SARS-CoV-2: Mutation_Count')
      if(is.null(outdir) == FALSE){
        ggsave(p,filename = paste(country, "_SNP_count.png", sep = ""),width = 12, height = 8, dpi=300, path = outdir)
      }
      if(is.null(outdir) == TRUE){
        print(p)
      }

    }
  }

}
