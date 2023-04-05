#' Reassign cells using knn
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
#'
#' @return A SingleCellExperiment with updated group assignments called 'knn'
#' @export
#'
#' @importFrom methods is
#' @importFrom class knn
#' @importFrom combinat combn
#' @import SingleCellExperiment
#'
#' @examples data(sce, snps)
#' sce <- high_conf_calls(sce)
#' sce <- add_snps(sce = sce, mat = snps, thresh = 0.8)
#' sce <- reassign(sce = sce, k = 10)
#'
reassign <- function(sce, 
                     k = 10, 
                     d = 10, 
                     train_cells = sce$train, 
                     predict_cells = sce$predict) {
    # Input checks
    stopifnot("'sce' must be of class SingleCellExperiment" = is(sce, "SingleCellExperiment"))
    stopifnot("k must be greater than or equal to two" = k > 1)
    stopifnot("k must be an integer" = k == round(k))

    # Singlet training data
    train <- counts(altExp(sce, "SNP"))[, train_cells == TRUE]
    labels <- sce$labels[train_cells == TRUE]
    labels <- droplevels(labels)
    colnames(train) <- labels

    # Simulated doublets
    combs <- expand.grid(levels(labels), levels(labels))
    combs <- combs[combs$Var1 != combs$Var2, ]
    p<-combn(unique(combs$Var1),2)
    combs_joined <- paste(p[1,],p[2,])

    all <- matrix(data=NA,
                  nrow=dim(altExp(sce,"SNP"))[1],
                  ncol=d*length(combs_joined))
    
    for (i in seq_along(combs_joined)) {
        l1 <- as.character(combs$Var1[i])
        l2 <- as.character(combs$Var2[i])

        d1 <- train[, labels == l1]
        d2 <- train[, labels == l2]

        s1 <- sample(seq_len(dim(d1)[2]), d , replace = TRUE)
        s2 <- sample(seq_len(dim(d2)[2]), d , replace = TRUE)

        doubs <- d1[, s1] == 1 | d2[, s2] == 1
        doubs <- doubs * 1
        doubs[doubs == 0] <- c(-1)

        all[,c(1+(i-1)*d):c(d+(i-1)*d)]<-doubs
        colnames(all)<-rep("Doublet",dim(all)[2])
    }
    train_all <- cbind(train, all)
    # prediction data
    pred <- as.data.frame(counts(altExp(sce, "SNP"))[, predict_cells == TRUE])
    # knn reclassification
    ID <- knn(t(train_all), t(pred), k = k, colnames(train_all))
    sce$knn <- as.character(sce$labels)
    sce$knn[predict_cells == TRUE] <- as.character(as.factor(ID))
    sce$knn <- as.factor(sce$knn)
    return(sce)
}
