data(sce,snps, package = "demuxSNP")

test_that("Only accepts SCE object", {
    seurat <- Seurat::as.Seurat(sce, data = NULL)
    expect_error(reassign(seurat))
})


test_that("Check reproducibility of knn", {
  sce<-add_snps(sce,snps)
  small_sce<-sce[,1:100]
  small_sce<-high_conf_calls(small_sce)
  set.seed(1)
  small_sce<-reassign(small_sce)
  expect_equal(sum(as.numeric(small_sce$knn)),356)
})