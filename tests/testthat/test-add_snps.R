test_that("Only accepts SCE object", {
    data(sce, snps, package = "demuxSNP")
    expect_error(add_snps(snps, snps))
})
