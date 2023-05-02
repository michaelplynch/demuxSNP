data(multiplexed_scrnaseq_sce, vartrix_consensus_snps, package = "demuxSNP")

test_that("Only accepts SCE object", {
    expect_error(add_snps(vartrix_consensus_snps, vartrix_consensus_snps))
})

test_that("Test reproducibility of SNP matrix filtering", {
    multiplexed_scrnaseq_sce<-add_snps(multiplexed_scrnaseq_sce, vartrix_consensus_snps)
    expect_equal(sum(counts(altExp(multiplexed_scrnaseq_sce,"SNP"))),c(-13423))
})
