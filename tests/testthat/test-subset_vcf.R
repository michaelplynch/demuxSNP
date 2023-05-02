data(multiplexed_scrnaseq_sce,commonvariants_1kgenomes_subset, package = "demuxSNP")

test_that("Only accept character genes", {
    top_genes <- seq_len(5)
    ensdb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
    expect_error(subset_vcf(commonvariants_1kgenomes_subset, top_genes, ensdb))
})

test_that("Check reproducibility in vcf subsetting", {
  top_genes <- common_genes(multiplexed_scrnaseq_sce)
  ensdb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
  vcf_sub<-subset_vcf(commonvariants_1kgenomes_subset, top_genes, ensdb)
  expect_equal(dim(vcf_sub)[1],2399)
})