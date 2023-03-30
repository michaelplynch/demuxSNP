data(sce, package = "demuxSNP")

test_that("Only accepts SCE object", {
    seurat <- Seurat::as.Seurat(sce, data = NULL)
    expect_error(high_conf_calls(seurat))
})

test_that("Track any changes in demuxmix wrapper", {
    sce<-sce[,1:100]
    sce<-high_conf_calls(sce)
    expect_equal(sum(as.numeric(sce$labels)),557)
})