#' Return a character vector of top n most frequent genes from a 
#' SingleCellExperiment object.
#'
#' @description Returns a character vector of the top n most frequently 
#' expressed genes from the counts of the SingleCellExperiment object.
#' Expression is based on having a count > 0 in a given cell.
#'
#' @param sce a SingleCellExperiment object
#' @param n number of genes to be returned. Defaults to n=100.
#'
#' @return character vector of n most frequently expressed genes.
#' @export
#'
#' @importFrom SingleCellExperiment counts
#' @importFrom Matrix rowSums
#' @importFrom methods is
#'
#' @examples data(multiplexed_scrnaseq_sce)
#' multiplexed_scrnaseq_sce <- common_genes(multiplexed_scrnaseq_sce)
#'
common_genes <- function(sce, n = 100) {
    # Input checks
    stopifnot("'sce' must be of class SingleCellExperiment" = is(sce, "SingleCellExperiment"))

    # calculating number of cells in which a given gene has count >=1 and sort
    counts_matrix <- counts(sce)
    gene_present <- counts_matrix > 0
    prop_cells <- rowSums(gene_present) / ncol(gene_present)
    sorted_top_genes <- sort(prop_cells, decreasing = TRUE)[seq_len(n)]
    top_genes <- names(sorted_top_genes)
    return(top_genes)
}
