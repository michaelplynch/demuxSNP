#' Subset common variants vcf file to only SNPs seen in most commonly expressed genes
#'
#' @param my_vcf path to SNPs vcf file
#' @param top_genes output from 'common_genes' function
#' @param ensdb object of class ensdb corresponding to organism, genome of data
#' @return common_snps.vcf vcf file containing subset of SNPs from common variants file seen in commonly expressed genes
#' @export
#' @import ensembldb
#' @import VariantAnnotation
#' @examples top_genes <- common_genes(sce)
#' ensdb<-EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
#' small_vcf <- subset_vcf(vcf, top_genes,ensdb)
#'
subset_vcf <- function(my_vcf, top_genes, ensdb) {
    SNP_ranges <- SummarizedExperiment::rowRanges(my_vcf)

    my_vcf_inbound <- my_vcf[BiocGenerics::end(SNP_ranges) <= GenomeInfoDb::seqlengths(SNP_ranges)[as.character(GenomeInfoDb::seqnames(SNP_ranges))]]
    SNP_ranges_inbound <- SummarizedExperiment::rowRanges(my_vcf_inbound)

    gns <- ensembldb::genes(ensdb)#EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86

    top_gene_ranges <- gns[gns$gene_name %in% top_genes]

    GenomeInfoDb::seqlengths(SNP_ranges_inbound) <- NA
    GenomeInfoDb::seqlengths(top_gene_ranges) <- NA

    top_genes_vcf <- my_vcf_inbound[IRanges::overlapsAny(SNP_ranges_inbound, top_gene_ranges, type = "within")]

    return(top_genes_vcf)
}
