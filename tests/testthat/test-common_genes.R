data(multiplexed_scrnaseq_sce, package = "demuxSNP")

test_that("Only accepts SCE object", {
    seurat <- Seurat::as.Seurat(multiplexed_scrnaseq_sce, data = NULL)
    expect_error(common_genes(seurat))
})

test_that("check consistency of common_genes calculation", {
  expect_equal(common_genes(multiplexed_scrnaseq_sce)[1:5],c("TPT1", "RPL13", "RPL28", "TMSB4X", "RPS27"))
})
