#' Add SNPs to sce
#'
#' @param sce Singnle Cell Experiment object
#' @param mat .mxt output from VartTrix
#' @param thresh threshold presence of cells
#'
#' @return
#' @export
#'
#' @examples
#'
add_snps<-function(sce,mat,thresh=0.8) {
  mat_obs<-mat[(rowSums(mat>0)/dim(mat)[2]) > thresh,]
  se<-SingleCellExperiment::SingleCellExperiment(list(counts=mat_obs))
  colnames(se)<-colnames(sce)
  rownames(se)<-rep("Snp",dim(mat_obs)[1])
  SingleCellExperiment::altExp(sce,"SNP")<-se
  return(sce)
}
