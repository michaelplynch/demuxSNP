#' Run existing cell hashing algorithms
#'
#' @param sce object of class SingleCellExperiment with HTO altExp assay
#'
#' @return Updated SingleCellExperiment object with logical vector indicating training data, data to be classified (complement of training data) and labels.
#' @export
#'
#' @importFrom methods is
#'
#' @examples sce <- consensus_calls(sce)
#'
consensus_calls <- function(sce) {
    ##Check input
    stopifnot("'sce' must be of class SingleCellExperiment"=is(sce,"SingleCellExperiment"))

    ##create training data
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
