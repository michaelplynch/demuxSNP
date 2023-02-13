#' Run existing cell hashing algorithms
#'
#' @param sce SingleCellExperiment object
#'
#' @return Updated SingleCellExperiment object
#' @export
#'
#' @examples
consensus_calls<-function(sce) {
# CiteFuse
result_sce<-CiteFuse::crossSampleDoublets(sce)

# Seurat
SingleCellExperiment::altExp(sce,"SNP")<-NULL
result_seurat<-Seurat::as.Seurat(sce)
#result_seurat$ident<-NULL
result_seurat<-Seurat::HTODemux(result_seurat)

# append results
sce$citefuse<-factor(result_sce$doubletClassify_between_label,levels=c("1","2","3","4","5","6","doublet/multiplet","negative"))
levels(sce$citefuse)<-c("Hashtag1","Hashtag2","Hashtag3","Hashtag4","Hashtag5","Hashtag6","Doublet","Negative")
sce$seurat<-factor(result_seurat$hash.ID,levels=c("Hashtag1","Hashtag2","Hashtag3","Hashtag4","Hashtag5","Hashtag6","Doublet","Negative"))
sce$intersect<-sce$citefuse==sce$seurat

sce$citefuse_class<-result_sce$doubletClassify_between_class
sce$seurat_class<-result_seurat$HTO_classification.global

sce$train<-sce$intersect==TRUE & sce$citefuse_class=="Singlet" & sce$seurat_class=="Singlet"
sce$predict<-sce$intersect==FALSE | sce$intersect==TRUE#& sce$citefuse_class=="Singlet" & sce$seurat_class=="Singlet"


return(sce)
}
