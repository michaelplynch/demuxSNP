#' Run existing cell hashing algorithms
#'
#' @param sce SingleCellExperiment object
#'
#' @return Updated SingleCellExperiment object
#' @export
#'
#' @examples sce <- consensus_calls(sce)
#'
consensus_calls <- function(sce) {
    rna <- BiocGenerics::colSums(SingleCellExperiment::counts(sce) > 0)
    hto <- as.matrix(SingleCellExperiment::counts(SingleCellExperiment::altExp(sce, "HTO")))

    dmm <- demuxmix::demuxmix(hto, rna = rna)
    demuxmix::pAcpt(dmm) <- 0.95
    classes <- demuxmix::dmmClassify(dmm)
    sce$train <- classes$Type == "singlet"
    sce$predict <- rep(TRUE, length(sce$train))
    sce$labels <- as.factor(classes$HTO)

    return(sce)
}
