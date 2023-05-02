data(multiplexed_scrnaseq_sce, package = "demuxSNP")

test_that("Only accepts SCE object", {
    seurat <- Seurat::as.Seurat(multiplexed_scrnaseq_sce, data = NULL)
    expect_error(high_conf_calls(seurat))
})

test_that("Track any changes in demuxmix wrapper", {
    sce<-multiplexed_scrnaseq_sce[,1:100]
    sce<-high_conf_calls(sce)
    expect_equal(sum(as.numeric(sce$labels)),557)
})