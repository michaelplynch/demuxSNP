# library(Matrix)
# library(SingleCellExperiment)
# #library(demuxSNP)
# data(multiplexed_scrnaseq_sce, vartrix_consensus_snps)
# sce<-multiplexed_scrnaseq_sce
# snps<-vartrix_consensus_snps
# 
# sce <- high_conf_calls(sce)
# sce <- add_snps(sce = sce, mat = snps, thresh = 0.8)
# 
# #reassign_balanced <- function(sce, k = 10, d = 10, train_cells = sce$train, predict_cells = sce$predict) {
# sce
# k = 10
# d = 10
# train_cells = sce$train
# predict_cells = sce$predict  
# 
# # Input checks
#   stopifnot("'sce' must be of class SingleCellExperiment" = is(sce, "SingleCellExperiment"))
#   stopifnot("k must be greater than or equal to two" = k > 1)
#   stopifnot("k must be an integer" = k == round(k))
#   
#   # rebalance training 
#   n <- min(table(sce$labels[train_cells]))
#   print(n)
#   if (n<50) {
#     n=50
#   }
#   print(n)
#   df<-data.frame(colData(sce[,train_cells]))
#   df$barcodes<-rownames(df)
#   df_balanced <- df %>% group_by(labels) %>% slice_sample(n=n)
#   
#   # Singlet training data
#   train <- counts(altExp(sce, "SNP"))[, colnames(sce) %in% df_balanced$barcodes]
#   train[train == -1] <- 0
#   labels <- sce$labels[colnames(sce) %in% df_balanced$barcodes]
#   labels <- droplevels(labels)
#   colnames(train) <- labels
#   combs <- expand.grid(levels(labels), levels(labels))
#   combs <- combs[combs$Var1 != combs$Var2, ]
#   p <- combn(unique(combs$Var1), 2)
#   combs_joined <- paste(p[1, ], p[2, ])
#   all <- matrix(data = NA, nrow = dim(altExp(sce, "SNP"))[1], 
#                 ncol = d * length(combs_joined))
#   for (i in seq_along(combs_joined)) {
#     l1 <- as.character(combs$Var1[i])
#     l2 <- as.character(combs$Var2[i])
#     d1 <- train[, labels == l1]
#     d2 <- train[, labels == l2]
#     s1 <- sample(seq_len(dim(d1)[2]), d, replace = TRUE)
#     s2 <- sample(seq_len(dim(d2)[2]), d, replace = TRUE)
#     doubs <- d1[, s1] == 1 | d2[, s2] == 1
#     doubs <- doubs * 1
#     doubs[doubs == 0] <- c(0)
#     all[, c(1 + (i - 1) * d):c(d + (i - 1) * d)] <- doubs
#     colnames(all) <- rep("Doublet", dim(all)[2])
#   }
#   train_all <- cbind(train, all)
#   pred <- as.data.frame(counts(altExp(sce, "SNP"))[, predict_cells == TRUE])
#   pred[pred == -1] <- 0
#   y = as.integer(as.factor(colnames(train_all)))
#   result_kern <- KernelKnn::KernelKnn(data = t(train_all), TEST_data = t(pred), 
#                                       y = y, k = k, Levels = seq_along(levels(as.factor(colnames(train_all)))), regression = F, method = "jaccard_coefficient")
#   ID <- apply(result_kern, 1, which.max)
#   z <- levels(as.factor(colnames(train_all)))
#   names(z)<-seq_along(z)
#   ID<-recode(as.character(ID), !!!z)
#   
#   sce$knn_jacc <- as.character(sce$labels)
#   sce$knn_jacc[predict_cells == TRUE] <- as.character(as.factor(ID))
#   sce$knn_jacc <- as.factor(sce$knn_jacc)
#   
#   print(dim(train))
#   print(dim(all))
#   
# #  return(sce)
# #}
# 
