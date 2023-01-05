common_genes<-function(sce) {
  counts_matrix<-SingleCellExperiment::counts(sce)
  gene_present<-counts_matrix>1
  prop_cells<-Matrix::rowSums(gene_present)/ncol(gene_present)
  sorted_top_genes<-base::sort(prop_cells,decreasing = TRUE)[1:100]
  top_genes<-base::names(sorted_top_genes)
  return(top_genes)
}
