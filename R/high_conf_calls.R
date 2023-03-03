#' Run existing cell hashing algorithms
#'
#' @param sce Object of class SingleCellExperiment with HTO altExp assay
#' @param assay Name of altExp for cell hashing counts to be retrieve from
#'
#' @return Updated SingleCellExperiment object with logical vector indicating training data, data to be classified (complement of training data) and labels.
#' @export
#'
#' @importFrom methods is
#'
#' @examples sce <- high_conf_calls(sce)
#'
high_conf_calls <- function(sce,assay="HTO") {
    ##Check input
    stopifnot("'sce' must be of class SingleCellExperiment"=is(sce,"SingleCellExperiment"))

    ##create training data
    rna <- BiocGenerics::colSums(SingleCellExperiment::counts(sce) > 0)
    hto <- as.matrix(SingleCellExperiment::counts(SingleCellExperiment::altExp(sce, assay)))

    #run demuxmix
    dmm <- demuxmix::demuxmix(hto, rna = rna)
    demuxmix::pAcpt(dmm) <- 0.95
    classes <- demuxmix::dmmClassify(dmm)

    #add labels back on to sce object
    sce$train <- classes$Type == "singlet"
    sce$predict <- rep(TRUE, length(sce$train))

    sce$labels <- classes$HTO
    sce$labels[classes$Type=="multiplet"]<-"multiplet"
    sce$labels<-as.factor(sce$labels)

    return(sce)
}
