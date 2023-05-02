data(multiplexed_scrnaseq_sce,vartrix_consensus_snps, package = "demuxSNP")

test_that("Only accepts SCE object", {
    seurat <- Seurat::as.Seurat(multiplexed_scrnaseq_sce, data = NULL)
    expect_error(reassign(seurat))
})

test_that("Check reproducibility of knn", {
  multiplexed_scrnaseq_sce<-add_snps(multiplexed_scrnaseq_sce, vartrix_consensus_snps)
  small_sce<-multiplexed_scrnaseq_sce[,1:100]
  small_sce<-high_conf_calls(small_sce)
  set.seed(1)
  small_sce<-reassign(small_sce)
  expect_equal(sum(as.numeric(small_sce$knn)),349)
})