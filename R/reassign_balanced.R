#' Reassign cells using balanced knn with jaccard distance
#'
#' @description k-nearest neighbour classification of cells.
#' Training data is intended to be labels of cells confidently called using 
#' cell hashing based methods and their corresponding SNPs.
#' Prediction data can be remaining cells but can also include the training 
#' data.
#' Doublets are simulated by randomly combining 'd' SNP profiles from each 
#' grouping combination.
#'
#' @param sce object of class SingleCellExperiment
#' @param k number of neighbours used in knn, defaults to 10
#' @param d number of doublets per group combination to simulate, defaults to 10
#' @param train_cells logical vector specifying which cells to use to train 
#' classifier
#' @param predict_cells logical vector specifying which cells to classify
#' @param n number of cells per group (otherwise will be calculated from data)
#' @param nmin min n per class (where available)
#' @return A SingleCellExperiment with updated group assignments called 'knn'
#' @export
#'
#' @importFrom methods is
#' @importFrom KernelKnn KernelKnn
#' @importFrom utils combn
#' @importFrom dplyr recode slice_sample group_by %>%
#' @import SingleCellExperiment
#'
#' @examples data(multiplexed_scrnaseq_sce, vartrix_consensus_snps)
#' multiplexed_scrnaseq_sce <- high_conf_calls(multiplexed_scrnaseq_sce)
#' multiplexed_scrnaseq_sce <- add_snps(sce = multiplexed_scrnaseq_sce, 
#' mat = vartrix_consensus_snps, 
#' thresh = 0.8)
#' multiplexed_scrnaseq_sce <- reassign_balanced(sce = multiplexed_scrnaseq_sce, k = 10)
#'
reassign_balanced <- function(sce, k = 10, d = 10, train_cells = sce$train, predict_cells = sce$predict, n = NULL, nmin = 50) {
  # Input checks
  stopifnot("'sce' must be of class SingleCellExperiment" = is(sce, "SingleCellExperiment"))
  stopifnot("k must be greater than or equal to two" = k > 1)
  stopifnot("k must be an integer" = k == round(k))
  
  # rebalance training
  if (is.null(n)) {
    n <- min(table(as.character(sce$labels[train_cells])))
  }
  print(n)
  if (n<50) {
    n=nmin
  }
  print(n)
  df<-data.frame(colData(sce[,train_cells]))
  df$barcodes<-rownames(df)
  df_balanced <- df %>% group_by(labels) %>% slice_sample(n=n)
  print(df_balanced$labels)
  # Singlet training data
  train <- counts(altExp(sce, "SNP"))[, colnames(sce) %in% df_balanced$barcodes]
  train[train == -1] <- 0
  labels <- sce$labels[colnames(sce) %in% df_balanced$barcodes]
  labels <- droplevels(labels)
  colnames(train) <- labels
  combs <- expand.grid(levels(labels), levels(labels))
  combs <- combs[combs$Var1 != combs$Var2, ]
  p <- combn(unique(combs$Var1), 2)
  combs_joined <- paste(p[1, ], p[2, ])
  all <- matrix(data = NA, nrow = dim(altExp(sce, "SNP"))[1], 
                ncol = d * length(combs_joined))
  for (i in seq_along(combs_joined)) {
    l1 <- as.character(combs$Var1[i])
    l2 <- as.character(combs$Var2[i])
    d1 <- train[, labels == l1]
    d2 <- train[, labels == l2]
    s1 <- sample(seq_len(dim(d1)[2]), d, replace = TRUE)
    s2 <- sample(seq_len(dim(d2)[2]), d, replace = TRUE)
    doubs <- d1[, s1] == 1 | d2[, s2] == 1
    doubs <- doubs * 1
    doubs[doubs == 0] <- c(0)
    all[, c(1 + (i - 1) * d):c(d + (i - 1) * d)] <- doubs
    colnames(all) <- rep("Doublet", dim(all)[2])
  }
  train_all <- cbind(train, all)
  pred <- as.data.frame(counts(altExp(sce, "SNP"))[, predict_cells == TRUE])
  pred[pred == -1] <- 0
  y = as.integer(as.factor(colnames(train_all)))
  result_kern <- KernelKnn::KernelKnn(data = t(train_all), TEST_data = t(pred), 
                                      y = y, k = k, Levels = seq_along(levels(as.factor(colnames(train_all)))), regression = F, method = "jaccard_coefficient")
  ID <- apply(result_kern, 1, which.max)
  z <- levels(as.factor(colnames(train_all)))
  names(z)<-seq_along(z)
  ID<-recode(as.character(ID), !!!z)
  
  sce$knn_balanced <- as.character(sce$labels)
  sce$knn_balanced[predict_cells == TRUE] <- as.character(as.factor(ID))
  sce$knn_balanced <- as.factor(sce$knn_balanced)
  
  print(dim(train))
  print(dim(all))
  
  return(sce)
}
