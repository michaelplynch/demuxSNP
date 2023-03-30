test_that("Only accepts SCE object", {
    data(sce, snps, package = "demuxSNP")
    expect_error(add_snps(snps, snps))
})

test_that("Test reprodicibility of SNP matrix filtering", {
    data(sce, snps, package = "demuxSNP")
    sce<-add_snps(sce, snps)
    expect_equal(sum(counts(altExp(sce,"SNP"))),c(-13423))
})
