#' Assays for mutation detection using different primers and probes
#'
#' These assays include the primer sequences and the detection ranges in which
#' mutations may occur.
#'
#' @docType data
#'
#' @usage data(assays)
#'
#' @format A dataframe with 10 rows and 7 columns.
#'
#' @keywords datasets
#'
#' @source This data is created by Zhanglab in Xiamen University.
#'
#' @examples
#' data(assays)
"assays"


#' A list of places in China
#'
#' The list is used for displacing some original cities' names with "China" in order
#' to make the downstream analysis easier.
#'
#' @docType data
#'
#' @usage data(chinalist)
#'
#' @format A dataframe with 31 rows and 1 column.
#'
#' @keywords datasets
#'
#' @source This data is created by Zhanglab in Xiamen University.
#'
#' @examples
#' data(chinalist)
"chinalist"

#' Mutation annotation results produced by "indelSNP" function
#'
#' A large matrix which could be used for downstream analysis like mutation
#' statistics description.
#'
#' @docType data
#'
#' @usage data(covid_annot)
#'
#' @format A large matrix.
#'
#' @keywords datasets
#'
#' @source \url{https://www.gisaid.org/}
#' @examples
#' data(covid_annot)
"covid_annot"

#' "GFF3" format gene position data for SARS-Cov-2
#'
#' This "GFF3" data is used for counting the mutations in each gene in virus sample.
#'
#' @docType data
#'
#' @usage data(gene_position)
#'
#' @format A dataframe with 26 rows and 10 columns.
#'
#' @keywords datasets
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/}
#'
#' @examples
#' data(gene_position)
"gene_position"

#' "GFF3" format annotation data for SARS-Cov-2
#'
#' This "GFF3" data is used for annotating the effects of mutations in virus sample.
#'
#' @docType data
#'
#' @usage data(gff3)
#'
#' @format A dataframe with 26 rows and 10 columns.
#'
#' @keywords datasets
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=2697049}
#'
#' @examples
#' data(gff3)
"gff3"

#' Mutation information derived from "nucmer" SNP analysis
#'
#' The "nucmer.snpss" variant file is obtained by processing the SARS-Cov-2 genomics sequence
#' from Gisaid website(complete, high coverage only, low coverage exclusion, Host=human, Virus name=hCoV-19)
#' by "seqkit" software and "nucmer SNP-calling" scripts. The total sample size is 5465.
#'
#' @docType data
#'
#' @usage data(nucmer)
#'
#' @format A dataframe with 43730 rows(mutation sites) and 14 columns.
#'
#' @keywords datasets
#'
#' @source \url{https://www.gisaid.org/}
#'
#' @examples
#' data(nucmer)
"nucmer"

#' Preprocessed "nucmer.snpss" file using "nucmerRMD" function
#'
#' A dataset contains some group information subtracted from the "nucmer"
#' object by "nucmerRMD" function in order to best describe the results.
#'
#' @docType data
#'
#' @usage data(nucmerr)
#'
#' @format A dataframe with 150000 rows(downsampled mutation sites) and 10
#' columns.
#'
#' @keywords datasets
#'
#' @source \url{https://www.gisaid.org/}
#'
#' @examples
#' data(nucmerr)
"nucmerr"

#' SARS-Cov-2 genomic reference sequence from NCBI
#'
#' This reference sequence is derived from "fasta" file, preprocessed by
#' "read.fasta" function(refseq<-read.fasta("NC_045512.2.fa",forceDNAtolower=FALSE)[[1]]). It is used for annotating mutations in virus samples.
#'
#' @docType data
#'
#' @usage data(refseq)
#'
#' @format "SeqFastadna" characters.
#'
#' @keywords datasets
#'
#' @source \url{https://pubmed.ncbi.nlm.nih.gov/32015508/}
#'
#' @examples
#' data(refseq)
"refseq"

