test_that("Only accepts SCE object", {
  seurat<-Seurat::as.Seurat(sce,data=NULL)
  expect_error(high_conf_calls(seurat))
})
