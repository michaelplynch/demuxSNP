#' Run demuxmix to determine high-confidence calls
#'
#' @param sce Object of class SingleCellExperiment with HTO (or similar) altExp
#'  assay
#' @param assay Name of altExp for cell hashing counts to be retrieved from
#'
#' @return Updated SingleCellExperiment object with logical vector indicating 
#' training data, data to be classified (all cells) and assigned labels for all 
#' cells.
#' @export
#'
#' @importFrom methods is
#' @importFrom MatrixGenerics colSums
#' @import demuxmix
#' @import SingleCellExperiment
#'
#' @examples data(multiplexed_scrnaseq_sce)
#' multiplexed_scrnaseq_sce <- high_conf_calls(multiplexed_scrnaseq_sce)
#'
high_conf_calls <- function(sce, assay = "HTO") {
    ## Check input
    stopifnot("'sce' must be of class SingleCellExperiment" = is(sce, "SingleCellExperiment"))

    ## create training data
    rna <- colSums(counts(sce) > 0)
    hto <- as.matrix(counts(altExp(sce, assay)))

    # run demuxmix
    dmm <- demuxmix(hto, rna = rna)
    pAcpt(dmm) <- 0.95
    classes <- dmmClassify(dmm)

    # add labels back on to sce object
    sce$train <- classes$Type == "singlet"
    sce$predict <- rep(TRUE, length(sce$train))

    sce$labels <- classes$HTO
    sce$labels[classes$Type == "multiplet"] <- "multiplet"
    sce$labels <- as.factor(sce$labels)

    return(sce)
}
