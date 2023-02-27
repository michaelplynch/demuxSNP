#' Subset common variants vcf file to only SNPs seen in most commonly expressed genes
#'
#' @param vcf path to SNPs vcf file
#' @param top_genes output from 'common_genes' function
#' @param ensdb object of class ensdb corresponding to organism, genome of data
#' @return common_snps.vcf vcf file containing subset of SNPs from common variants file seen in commonly expressed genes
#' @export
#'
#' @import ensembldb
#' @import VariantAnnotation
#' @importFrom("methods", "is")
#'
#' @examples top_genes <- common_genes(sce)
#' ensdb<-EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
#' small_vcf <- subset_vcf(vcf, top_genes,ensdb)
#'
subset_vcf <- function(vcf, top_genes, ensdb) {
    #Input checks
    stopifnot("'vcf' must be of class CollapsedVCF"=is(vcf,"CollapsedVCF"))
    stopifnot("'top_genes' must be of class character"=is(top_genes,"character"))
    stopifnot("'ensdb' must be of class EnsDb"=is(ensdb,"EnsDb"))
    stopifnot("Are your SNPs and ensdb object from the same genome build?"=seqinfo(vcf)@genome[1]==seqinfo(ensdb)@genome[1])

    #calculating row ranges of SNPs, genes, then subsetting vcf by gene ranges
    SNP_ranges <- SummarizedExperiment::rowRanges(vcf)

    vcf_inbound <- vcf[BiocGenerics::end(SNP_ranges) <= GenomeInfoDb::seqlengths(SNP_ranges)[as.character(GenomeInfoDb::seqnames(SNP_ranges))]]
    SNP_ranges_inbound <- SummarizedExperiment::rowRanges(vcf_inbound)

    gns <- ensembldb::genes(ensdb)

    top_gene_ranges <- gns[gns$gene_name %in% top_genes]

    GenomeInfoDb::seqlengths(SNP_ranges_inbound) <- NA
    GenomeInfoDb::seqlengths(top_gene_ranges) <- NA

    top_genes_vcf <- vcf_inbound[IRanges::overlapsAny(SNP_ranges_inbound, top_gene_ranges, type = "within")]

    return(top_genes_vcf)
}
