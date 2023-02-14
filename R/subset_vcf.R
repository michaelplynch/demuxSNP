#' Subset common variants vcf file to only SNPs seen in commonly expressed genes
#'
#' @param my_vcf path to common variants vcf file
#' @param sce single cell experiment object
#' @param top_genes from previous function
#'
#' @return common_snps.vcf vcf file containing subset of SNPs from common variants file seen in commonly expressed genes
#' @export
#'
#'
#'
subset_vcf<-function(my_vcf,sce,top_genes) {
  #my_vcf<-VariantAnnotation::readVcf('full.vcf',"GRCh38")
  SNP_ranges<-SummarizedExperiment::rowRanges(my_vcf)

  my_vcf_inbound<-my_vcf[BiocGenerics::end(SNP_ranges)<=GenomeInfoDb::seqlengths(SNP_ranges)[as.character(GenomeInfoDb::seqnames(SNP_ranges))]]
  SNP_ranges_inbound<-SummarizedExperiment::rowRanges(my_vcf_inbound)

  gns <- ensembldb::genes(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)

  top_gene_ranges<-gns[mcols(gns)$gene_name %in% top_genes]

  GenomeInfoDb::seqlengths(SNP_ranges_inbound)<-NA
  GenomeInfoDb::seqlengths(top_gene_ranges)<-NA
  logi<-SNP_ranges_inbound %within% top_gene_ranges
  top_genes_vcf<-my_vcf_inbound[SNP_ranges_inbound %within% top_gene_ranges]
  top_genes_vcf
  #VariantAnnotation::writeVcf(top_genes_vcf,'reduced_SNPs.vcf')
  return(top_genes_vcf)
}

