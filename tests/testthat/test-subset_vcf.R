test_that("Only accept character genes", {
    data(vcf, package="demuxSNP")
    top_genes<-seq_len(5)
    ensdb<-EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
    expect_error(subset_vcf(vcf,top_genes,ensdb))
})
