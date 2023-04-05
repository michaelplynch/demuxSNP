#' SingleCellExperiment object containing multiplexed RNA and HTO data from
#' six biological smamples
#'
#' Example SingleCellExperiment object containing demultiplexed scRNAseq data 
#' from six donors, used throughout and built upon in demuxSNP workflow.
#'
#' @return ## `sce`
#' An object of class SingleCellExperiment
#' @usage  data(sce)
"sce"
#'
#' Sample VarTrix output
#'
#' A sample output from VarTrix corresponding to the sce SingleCellExperiment 
#' objec for a subset of SNPs located in well observed genes.
#'
#' @return ## `snps`
#' An object of class matrix
#' @usage data(snps)
"snps"
#'
#' Sample vcf file
#'
#' VCF file containing SNPs from a subset of the 1k Genomes common variants 
#' HG38 genome build.
#'
#' @return ## `vcf`
#' An object of class CollapsedVcf
#' @source https://cellsnp-lite.readthedocs.io/en/latest/snp_list.html
#' @usage data(vcf)
"vcf"