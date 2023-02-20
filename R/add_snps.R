#' Add SNPs to SingleCellExperiment object
#'
#' @param sce SingleCellExperiment object
#' @param mat .mtx output from VartTrix
#' @param thresh threshold presence of SNP, defaults to 0.8
#'
#' @return Updated SingleCellExperiment object
#' @export
#'
#' @examples sce<-add_snps(sce=sce,mat=snps,thresh=0.8)
#'
add_snps<-function(sce,mat,thresh=0.8) {
    mat_obs<-mat[(rowSums(mat>0)/dim(mat)[2]) > thresh,]

    mat_obs[mat_obs==0]<-c(0)
    mat_obs[mat_obs==1]<-c(-1)
    mat_obs[mat_obs==2]<-1
    mat_obs[mat_obs==3]<-1

    se<-SingleCellExperiment::SingleCellExperiment(list(counts=mat_obs))
    colnames(se)<-colnames(sce)
    rownames(se)<-rep("Snp",dim(mat_obs)[1])
    SingleCellExperiment::altExp(sce,"SNP")<-se
    return(sce)
}
