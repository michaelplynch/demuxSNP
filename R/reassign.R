#' Reassign cells using knn
#'
#' @param sce Single Cell Experiment object
#' @param k number of neighbours in knn, default to 10
#' @param seed set seed for reproducibility, default to 1
#' @param train_cells logical vector specifying which cells to use to train classifier
#' @param predict_cells logical vector specifying which cells to classify
#'
#' @return A SingleCellExperiment with updated group assignments in "knn" metadata
#' @export
#'
#' @examples
#' sce<-consensus_calls(sce)
#' sce<-add_snps(sce=sce,mat=snps,thresh=0.8)
#' sce<-reassign(sce=sce,k=10)
reassign<-function(sce,k=10,seed=1,train_cells=sce$train,predict_cells=sce$predict) {

  #Singlet training data
  train<-SingleCellExperiment::counts(SingleCellExperiment::altExp(sce,"SNP"))[,sce$train==TRUE]
  labels<-sce$labels[train_cells==TRUE]
  labels<-droplevels(labels)

  colnames(train)<-labels

  #Simulated doublets
  combs<-expand.grid(levels(labels),levels(labels))
  combs<-combs[combs$Var1!=combs$Var2,]
  combs_joined<-paste0(combs$Var1,combs$Var2)

  all<-c()
  for (i in seq_len(length(combs_joined))) {

    l1<-as.character(combs$Var1[i])
    l2<-as.character(combs$Var2[i])

    d1<-train[,labels==l1]
    d2<-train[,labels==l2]

    s1<-sample(1:dim(d1)[2],10,replace=TRUE)
    s2<-sample(1:dim(d2)[2],10,replace=TRUE)

    doubs<-d1[,s1]==1 | d2[,s2]==1
    doubs<-doubs*1
    doubs[doubs==0]<-c(-1)
    colnames(doubs)<-rep("Doublet",length(colnames(doubs)))#rep(paste0(l1,l2),length(colnames(doubs)))

    all<-cbind(all,doubs)
  }

  dim(all)

  train_all<-cbind(train,all)

  #prediction data
  pred<-as.data.frame(SingleCellExperiment::counts(SingleCellExperiment::altExp(sce,"SNP"))[,sce$predict==TRUE])

  #set.seed(seed)
  ID<-class::knn(t(train_all),t(pred),k=k,colnames(train_all))

  #sce$knn[1:length(colnames(sce))]<-"unknown"
  sce$knn<-as.factor(ID)
  return(sce)
}
