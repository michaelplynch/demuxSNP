data(sce,vcf, package = "demuxSNP")

test_that("Only accept character genes", {
    top_genes <- seq_len(5)
    ensdb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
    expect_error(subset_vcf(vcf, top_genes, ensdb))
})

test_that("Check reproducibility in vcf subsetting", {
  top_genes <- common_genes(sce)
  ensdb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
  vcf_sub<-subset_vcf(vcf, top_genes, ensdb)
  expect_equal(dim(vcf_sub)[1],2399)
})