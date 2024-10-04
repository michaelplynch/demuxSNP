#' Reassign cells based on SNPs
#'
#' @param sce SingleCellExperiment object
#' @param train_cells logical, cells to be used for training
#' @param predict_cells logical, cells to be used for prediction
#' @param labels provisional cell labels
#' @param min_cells minimum coverage (number of cells with read at SNP location) for SNP to be used for classification.
#' @param key unique key in naming of singlet groups used with grep to remove doublet/negative/uncertain labels
#' @return character vector containing reassignments
#' @export
#' @importFrom stats aggregate
#' @examples
#' data(multiplexed_scrnaseq_sce, vartrix_consensus_snps)
#' multiplexed_scrnaseq_sce <- high_conf_calls(multiplexed_scrnaseq_sce)
#' multiplexed_scrnaseq_sce <- add_snps(sce = multiplexed_scrnaseq_sce, 
#' mat = vartrix_consensus_snps, 
#' thresh = 0.8)
#' multiplexed_scrnaseq_sce<-reassign_centroid(multiplexed_scrnaseq_sce)
reassign_centroid<-function(sce,train_cells=sce$train,predict_cells=sce$predict,labels=sce$labels,min_cells=30,key="Hashtag") {
  
  snps<-counts(altExp(sce,"SNPcons"))
  # find cluster multivariate modes
  agg<-cluster_modes(snps_cons = snps, training_labels = labels, min_cells = min_cells)
  #print(head(agg))
  print(dim(agg))
  print(sum(agg))
  # simulate doublets
  train_mat<-sim_doubs(centroids=agg[,grep(key,colnames(agg))])
  
  # create distance matrix
  jw<-jaccard_weighted(snps_predict=snps ,snps_train=train_mat)
  
  # perform knn
  nearest_ns<-knearestneighbours(jw,k=1,train_cells,predict_cells,train_labels = colnames(train_mat))
  nearest_ns[!predict_cells]<-as.character(labels[!predict_cells])
  #nearest_ns<-gsub('Hashtag','K',nearest_ns) ## come back to this and make it generalisable
  return(nearest_ns)
}


sim_doubs<-function(centroids) {
  agg_doub<-centroids
  labels <- colnames(agg_doub)
  p <- combn(unique(labels), 2)
  combs_joined <- paste(p[1, ], p[2, ])
  all <- matrix(data = NA, nrow = dim(centroids)[1],
                ncol = length(combs_joined))
  for (i in seq_along(combs_joined)) {
    l1 <- as.character(p[1,i])
    l2 <- as.character(p[2,i])
    d1 <- agg_doub[, labels == l1]
    d2 <- agg_doub[, labels == l2]
    
    d1[d2==0]<-0 # should help balance doublet profile being driven by singlet 
    d2[d1==0]<-0 # group with most SNPs. Takes the union of present SNPs
    doubs<-matrix(1,nrow = length(d1),ncol=1)
    doubs[d1 ==0 | d2 ==0]<-0
    #doubs[d1 %in% c(1) & d2 %in% c(1)]<-1
    doubs[d1 %in% c(2,3) | d2 %in% c(2,3)]<-3
    all[, i] <- doubs
  }
  #colnames(all) <- paste0("Doublet",p[1,],p[2,])
  colnames(all)<-rep('Doublet',ncol(all))
  train_all <- cbind(agg_doub, all)
  print(colSums(train_all))
  return(train_all)
  
}

jaccard_weighted<-function(snps_predict,snps_train) {
  mat<-t(snps_predict[])
  mat[mat==0]<-0
  mat[mat==1]<-0
  mat[mat==2]<-1
  mat[mat==3]<-1
  
  mat_bool<-t(snps_predict!=0)
  mat_bin <- mat_bool*1
  
  mat_train<-t(snps_train[])
  mat_train[mat_train==0]<-0
  mat_train[mat_train==1]<-0
  mat_train[mat_train==2]<-1
  mat_train[mat_train==3]<-1
  
  mat_bool_train<-t(snps_train!=0)
  mat_bin_train <- mat_bool_train*1
  
  #mat_bin<-matrix(1,ncol=ncol(mat_bin),nrow=nrow(mat_bin))
  #mat_bin_train<-matrix(1,ncol=ncol(mat_bin_train),nrow=nrow(mat_bin_train))
  a <- mat %*% t(mat_bin_train*mat_train)
  b <- (mat*mat_bin) %*% ((1 - t(mat_train))*t(mat_bin_train))
  c <- ((1 - mat)*mat_bin) %*% t(mat_train*mat_bin_train)
  d <- ((1 - mat)*mat_bin) %*% t((1 - mat_train)*mat_bin_train)
  
  #print(a[1:5,1:6])
  #print(b[1:5,1:6])
  #print(c[1:5,1:6])
  #print(d[1:5,1:6])
  
  j<-1-((a)/(a+b+c))
  return(j)
}

modeValue <- function(v) {
  if (length(v)==1) 
    return(v)
  #construct the frequency table
  tbl = table(v) 
  #get the most frequently occurring value 
  myMode = names(tbl)[which.max(tbl)]
  return (myMode=myMode)
}

knearestneighbours <- function(dist_mat,k,train_cells,predict_cells,train_labels) {
  
  trainlabs<-train_labels
  dist_adj<-dist_mat
  pred_labs<-c()
  for (ii in seq_along(predict_cells)) {
    if (ii%%1500 == 0) {
      print(ii)
    }
    nearestpoints<-sort(dist_adj[ii,],index.return=TRUE,na.last=T)$ix[1:k]
    nearestlabs<-trainlabs[nearestpoints]
    pred_labs[ii]<-modeValue(nearestlabs)
  }
  return(pred_labs)
}

Mode <- function(x,min) {
  ux <- unique(x[x!=0])
  mod<-ux[which.max(tabulate(match(x, ux)))]
  if (sum(x %in% mod)<min){
    return(0)
  }
  else {
    return(mod)}
}

cluster_modes<-function(snps_cons,training_labels,min_cells) {
  snps<-snps_cons
  snps_tf<-t(snps)
  agg<-t(aggregate(snps_tf,by=list(training_labels), Mode,min=min_cells))
  agg[is.na(agg)]<-0
  colnames(agg)<-agg[1,]
  agg<-agg[-c(1),]
  agg<-agg[,grep("Hashtag",colnames(agg))]
  mode(agg)<-'numeric'
  return(agg)
}