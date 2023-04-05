#' Add SNPs to SingleCellExperiment object
#'
#' @param sce object of class SingleCellExperiment
#' @param mat object of class matrix, output from VarTrix in 'consensus' mode 
#' (default)
#' @param thresh threshold presence of SNP, defaults to 0.8
#'
#' @return Updated SingleCellExperiment object with snps in altExp slot
#' @export
#'
#' @importFrom methods is
#' @import SingleCellExperiment
#'
#' @examples data(sce, snps)
#' sce <- add_snps(sce = sce, mat = snps, thresh = 0.8)
#'
add_snps <- function(sce, mat, thresh = 0.8) {
    # Input checks
    stopifnot("'sce' must be of class SingleCellExperiment" = is(sce, "SingleCellExperiment"))
    stopifnot("thresh must be between 0 and 1" = thresh < 1 & thresh > 0)
    stopifnot("SingleCellExperiment and snps matrix contain unequal number of cells" = dim(counts(sce))[2] == dim(mat)[2])
    stopifnot("Did you run VarTrix in the default 'consensus' mode?" = identical(mat, round(mat)))

    # Add snps to sce
    mat_obs <- mat[(rowSums(mat > 0) / dim(mat)[2]) > thresh, ]

    mat_obs[mat_obs == 0] <- c(0)
    mat_obs[mat_obs == 1] <- c(-1)
    mat_obs[mat_obs == 2] <- 1
    mat_obs[mat_obs == 3] <- 1

    se <- SingleCellExperiment(list(counts = mat_obs))
    colnames(se) <- colnames(sce)
    rownames(se) <- rep("Snp", dim(mat_obs)[1])
    altExp(sce, "SNP") <- se
    return(sce)
}
