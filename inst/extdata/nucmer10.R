#' The "nucmer10.snps" variant file is obtained by processing the SARS-Cov-2 genomics sequence
#' from Gisaid website(complete, high coverage only, low coverage exclusion, Host=human, Virus name=hCoV-19)
#' by "seqkit" software and "nucmer SNP-calling" scripts. The total sample size is 5465.
#'
#' @docType extdata
#'
#' @keywords datasets
#'
#' @source \url{https://www.gisaid.org/}
#' @details The full analysis for reproducing the "nucmer10.snps" file:
#' #software under linux: seqkit , nucmer
#'
#' #download the SARS-Cov-2 genomics sequence(*.fasta) from Gisaid website( comploe, high coverage only, low coverage exclusion, Host=human, Virus name=hCoV-19)
#' #download  NC_045512.2.fa from NCBI
#'
#' cd /mnt/e/2019_nCov/20200508/ #
#' input=gisaid_hcov-19_2020_06_14_04.fasta #your fasta file name.
#' ## clear the data
#' seqkit grep  -s -p - $input -v > Gisaid_clear.fasta  # remove the data with '-'
#' #Run nucmer to obtain variant file
#' ref=NC_045512.2.fa # The reference SARS-CoV-2 Wuhan Genome
#' ## Covert the DOS/window file format to UNIX format for both ref and input files
#' sed 's/^M$//' Gisaid_clear.fasta > Gisaid_clear_format.fasta
#' sed 's/^M$//' NC_045512.2.fa > ref.fa
#' ### remove fasta sequence with duplicated ID
#' awk '/^>/{f=!d[$1];d[$1]=1}f' Gisaid_clear_format.fasta > Gisaid_RMD.fasta
#' ### calculate total sample
#' grep -c '^>' Gisaid_RMD.fasta
#' ### downsampling fasta seq:
#' seqkit sample --proportion 0.15 Gisaid_RMD.fasta > Gisaid_RMD_15.fasta
#' ### Snap-calling
#' nucmer --forward -p nucmer ref.fa Gisaid_RMD.fasta
#' show-coords -r -c -l nucmer.delta > nucmer.coords
#' show-snps nucmer.delta -T -l > nucmer.snps

#' @examples
#' options(stringsAsFactors = F)
#' dir<- system.file("extdata", package = "CovidMutations")
#' list.files(dir)
#' nucmer<- read.delim(unzip(dir),as.is=TRUE,skip=4,header=FALSE)
"nucmer10.snps"

