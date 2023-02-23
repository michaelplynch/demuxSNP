#' Return a list of the top 100 most commonly expressed genes from a SingleCellExperiment object.
#'
#' @param sce A SingleCellExperiment object
#'
#' @return names of top 100 most commonly expressed genes.
#' @export
#' @importFrom SingleCellExperiment counts
#' @importFrom Matrix rowSums
#' @examples
#' sce <- common_genes(sce)
#'
common_genes <- function(sce) {
    counts_matrix <- counts(sce)
    gene_present <- counts_matrix > 0
    prop_cells <- Matrix::rowSums(gene_present) / ncol(gene_present)
    sorted_top_genes <- sort(prop_cells, decreasing = TRUE)[seq_len(100)]
    top_genes <- names(sorted_top_genes)
    return(top_genes)
}
