#' SingleCellExperiment object containing multiplexed RNA and HTO data from
#' six biological smamples
#'
#' Example SingleCellExperiment object containing demultiplexed scRNAseq data 
#' from six donors, used throughout and built upon in demuxSNP workflow.
#'
#' @return ## `multiplexed_scrnaseq_sce`
#' An object of class SingleCellExperiment
#' @usage  data(multiplexed_scrnaseq_sce)
"multiplexed_scrnaseq_sce"
#'
#' Sample VarTrix output
#'
#' A sample output from VarTrix corresponding to the sce SingleCellExperiment 
#' objec for a subset of SNPs located in well observed genes.
#'
#' @return ## `vartrix_consensus_snps`
#' An object of class matrix
#' @usage data(vartrix_consensus_snps)
"vartrix_consensus_snps"
#'
#' Sample vcf file
#'
#' VCF file containing SNPs from a subset of the 1k Genomes common variants 
#' HG38 genome build.
#'
#' @return ## `commonvariants_1kgenomes_subset`
#' An object of class CollapsedVcf
#' @source https://cellsnp-lite.readthedocs.io/en/latest/snp_list.html
#' @usage data(commonvariants_1kgenomes_subset)
"commonvariants_1kgenomes_subset"