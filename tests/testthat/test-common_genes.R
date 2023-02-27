test_that("Only accepts SCE object", {
  seurat<-Seurat::as.Seurat(sce,data=NULL)
  expect_error(top_genes(seurat))
})
