#' Reassign cells using knn
#'
#' @param sce Single Cell Experiment object
#' @param k number of neighbours in knn, default to 10
#' @param seed set seed for reproducibility, default to 1
#'
#' @return
#' @export
#'
#' @examples
reassign<-function(sce,k=10,seed=1,train=sce$train,predict=sce$predict) {

  train<-as.data.frame(t(SingleCellExperiment::counts(SingleCellExperiment::altExp(sce,"SNP"))[,sce$train==TRUE]))
  pred<-as.data.frame(t(SingleCellExperiment::counts(SingleCellExperiment::altExp(sce,"SNP"))[,sce$predict==TRUE]))



  set.seed(seed)
  ID<-class::knn(train,pred,k=10,sce$citefuse[sce$train==TRUE])

  #sce$knn[1:length(colnames(sce))]<-"unknown"
  sce$knn<-as.factor(ID)
  return(sce)
}
