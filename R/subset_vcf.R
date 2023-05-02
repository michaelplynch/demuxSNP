#' Subset common variants vcf file to only SNPs seen in 'top_genes'
#'
#' @param vcf object of class CollapsedVCF
#' @param top_genes output from 'common_genes' function, alternatively character vector containing custom gene names.
#' @param ensdb object of class EnsDb corresponding to organism, genome of data
#'
#' @return object of class CollapsedVCF containing subset of SNPs from supplied vcf seen in commonly expressed genes
#' @export
#'
#' @import ensembldb
#' @import VariantAnnotation
#' @importFrom methods is
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom BiocGenerics end
#' @importFrom IRanges overlapsAny
#' @import GenomeInfoDb
#'
#' @examples data(multiplexed_scrnaseq_sce, commonvariants_1kgenomes_subset)
#' top_genes <- common_genes(multiplexed_scrnaseq_sce)
#' ensdb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
#' small_vcf <- subset_vcf(commonvariants_1kgenomes_subset, top_genes, ensdb)
#'
subset_vcf <- function(vcf, top_genes, ensdb) {
    # Input checks
    stopifnot("'vcf' must be of class CollapsedVCF" = is(vcf, "CollapsedVCF"))
    stopifnot("'top_genes' must be of class character" = is(top_genes, "character"))
    stopifnot("'ensdb' must be of class EnsDb" = is(ensdb, "EnsDb"))
    stopifnot("Are your SNPs and ensdb object from the same genome build?" = genome(vcf)[1] == genome(ensdb)[1])

    # calculating row ranges of SNPs, genes, then subsetting vcf by gene ranges
    SNP_ranges <- rowRanges(vcf)

    vcf_inbound <- vcf[end(SNP_ranges) <= seqlengths(SNP_ranges)[as.character(seqnames(SNP_ranges))]]
    SNP_ranges_inbound <- rowRanges(vcf_inbound)

    gns <- ensembldb::genes(ensdb)

    top_gene_ranges <- gns[gns$gene_name %in% top_genes]

    seqlengths(SNP_ranges_inbound) <- NA
    seqlengths(top_gene_ranges) <- NA

    top_genes_vcf <- vcf_inbound[overlapsAny(SNP_ranges_inbound, top_gene_ranges, type = "within")]

    return(top_genes_vcf)
}
