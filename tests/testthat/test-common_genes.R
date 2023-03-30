test_that("Only accepts SCE object", {
    data(sce, package = "demuxSNP")
    seurat <- Seurat::as.Seurat(sce, data = NULL)
    expect_error(common_genes(seurat))
})

test_that("check consistency of common_genes calculation", {
  data(sce, package = "demuxSNP")
  expect_equal(common_genes(sce)[1:5],c("TPT1", "RPL13", "RPL28", "TMSB4X", "RPS27"))
})
