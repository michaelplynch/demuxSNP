test_that("Only accepts SCE object", {
    data(sce, package="demuxSNP")
    seurat<-Seurat::as.Seurat(sce,data=NULL)
    expect_error(reassign(seurat))
})
