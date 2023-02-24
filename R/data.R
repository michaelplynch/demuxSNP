#' SingleCellExperiment object containing RNA and HTO data
#'
#' Example SingleCellExperiment object used throughout and built upon in SNPcheck workflow.
#'
#' @format ## `sce`
#' A data frame with 7,240 rows and 60 columns:
#' @source n/a for now
#' @examples sce <- data(sce)
"sce"
#'
#' Sample VarTrix output
#'
#' A sample output from VarTrix corresponding to the sce SingleCellExperiment object.
#'
#' @format ## `snps`
#' A matrix with 2,000 cells and 2k odd SNPs
#' @source n/a for now
#' @examples dim(snps)
"snps"
#'
#' Sample vcf file
#'
#' A subset of the 1k Genomes common variants HG38 vcf file.
#'
#' @format ## `vcf`
#' A vcf file containing x variants.
#' @source <https://www.who.int/teams/global-tuberculosis-programme/data>
"vcf"
