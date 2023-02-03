#' Return a list of the top 100 most commonly expressed genes from a SingleCellExperiment object.
#'
#' @param sce A SingleCellExperiment object
#'
#' @return names of top 100 most commonly expressed genes.
#' @export
#'
#'
#'
common_genes<-function(sce) {
  counts_matrix<-SingleCellExperiment::counts(sce)
  gene_present<-counts_matrix>1
  prop_cells<-Matrix::rowSums(gene_present)/ncol(gene_present)
  sorted_top_genes<-base::sort(prop_cells,decreasing = TRUE)[1:100]
  top_genes<-base::names(sorted_top_genes)
  return(top_genes)
}
