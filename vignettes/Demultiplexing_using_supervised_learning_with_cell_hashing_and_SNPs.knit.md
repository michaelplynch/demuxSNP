---
title: "Demultiplexing_using_supervised_learning_with_cell_hashing_and_SNPs"
author:
  - name: Michael Lynch, University of Limerick
  - name: Aedin Culhane, University of Limerick
output:
  BiocStyle::html_document:
    toc_float: true
bibliography: references.bib
vignette: |
  %\VignetteIndexEntry{dim reduction with corral}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
  
---







```r

library(SNPcheck)
library(ComplexHeatmap)
library(viridisLite)
library(Seurat)
library(ggpubr)
library(dittoSeq)
library(utils)
```


```r
colors <- structure(viridis(n = 3), names = c("-1", "0", "1"))
```

# Introduction

Demultiplexing in scRNAseq involves assigning cells back to their original sample, where cells from different donors, treatment types or physiological locations are sequenced together.
This allows for a significant reduction in cost of sequencing as well as facilitating removal of technical artifacts such as batch effects and multi sample doublets.
To address this need, a large number of methods have been proposed. 
However, a universally robust algorithms remains elusive.
Below, we introduce some existing methods and highlight the novel features of our approach and its advantages to the user.

## Existing methods


### Cell Hashing

Cells from each group are labelled with a distinct tag (HTO or LMO) which is sequenced to give a counts matrix.
Due to non-specific binding, these counts generally form a bimodal distribution.
Such methods are generally computationally efficient.
Their performance, however, is highly dependent on the tagging quality.

Binary vs probabilistic methods

To extend this, recent methods now attempt to probabilistically model this process, allowing users to define a cut-off threshold for the assignment confidence.
While these methods give the user greater flexibility in determinng which cells to keep, the question remains as to what an appropriate value for a sensible cut off is.
This also results in the removal of potentially high quality, valuable cells from downstream analysis.

@boggy_bff_2022

@stoeckius_cell_2018

@kim_citefuse_2020



### SNPs

The second class of methods exploits natural genetic variation between cells and so can only be used where the groups are genetically distinct.
Demuxlet @kang_multiplexed_2018 -high accuracy but requires genotype information
Souporcell @heaton_souporcell_2020 and Vireo @huang_vireo_2019
The performance of these methods is reduced by the presence of ambient RNA and unequal donor contributions.

Demuxlet remains the standard used to benchmark other methods.

## Motivation

Here, we aim to leverage information from both modalities to optimise classification and minimise waste.
We note that as with other SNP based methods, genetic difference between sampe groups is a prerequisite, as well as sufficient sequencing depth.

Novel features:

* Uses both cell hashing and SNP data. 
Current methods are limited to using either the cell hashing counts or SNP calls.
By using a supervised approach, performance increases in classification, particularly when cell contributions are unevenly distributed.
* Selects SNPs based on gene expression to reduce noise and computational cost.
This reduces artifacts caused by cell-type and speeds things up.



## Installation


```r

devtools::install_github("michaelplynch/SNPcheck", build_vignettes = TRUE)

browseVignettes(package="SNPcheck")

```

## Quick Usage


```r
# subset common variants file:
top_genes<-common_genes(sce)
small_vcf<-subset_vcf(sce,vcf)

# create training (high confidence) data
sce<-consensus_calls(sce)

## Reassignment
sce<-add_snps(sce,mat)
sce<-reassign(sce)

```

# Exploratory analysis

We load three data objects. 
A SingleCellExperiment object containing RNA and HTO counts, a vcf file containing SNPs and a matrix containing SNP information for each cell (we will show you how to generate this SNPs matrix using VarTrix outside of R).


```r
data(sce,vcf,snps,package="SNPcheck")
```

The HTO or LMO distribution is usually bimodal, with a signal (high counts) and background distribution (low counts) caused by non-specific binding.
Ideally, these distributions would be clearly separated with no overlap, but in practice, this is not always the  case.
in our example data, we see that the signal and noise overlap to varying extents in each Hashtag.

<img src="C:/Users/michael.lynch/AppData/Local/Temp/RtmpMFq7J3/preview-46a04eb720d8.dir/Demultiplexing_using_supervised_learning_with_cell_hashing_and_SNPs_files/figure-html/unnamed-chunk-7-1.png" width="100%" />

As an example, we will run HTODemux on the data.


```r
seurat <- as.Seurat(sce,data=NULL)
seurat <- HTODemux(seurat)
#> Warning in PseudobulkExpression(object = object, pb.method = "average", :
#> Exponentiation yielded infinite values. `data` may not be log-normed.
seurat$hash.ID <- factor(as.character(seurat$hash.ID))
sce$seurat <- seurat$hash.ID

sce$seurat <- seurat$hash.ID

table(sce$seurat)
#> 
#>  Doublet Hashtag1 Hashtag2 Hashtag3 Hashtag4 Hashtag5 Hashtag6 Negative 
#>      328       15       23      298      189       12      500      635
```

Here, we see an unusually large number of cells being called as "Negative".

Additionally, the library size of the "Negative" group looks similar to that of other groups.



```r
seurat$libsize <- colSums(GetAssayData(seurat,slot="counts",assay="RNA"))
dittoPlot(seurat, "libsize", group.by = "hash.ID")
```

<img src="C:/Users/michael.lynch/AppData/Local/Temp/RtmpMFq7J3/preview-46a04eb720d8.dir/Demultiplexing_using_supervised_learning_with_cell_hashing_and_SNPs_files/figure-html/unnamed-chunk-9-1.png" width="100%" />

For the remainder of this vignette we'll outline our method of checking whether or not cells have been called correctly and how to assign them to their appropriate group!

# Preprocessing

Common variants files, for example from the 1000 Genomes Project, can contain over 7 million SNPs.
To reduce computational cost and cell-type effects, we subset our SNPs list to those located within genes expressed across most cells in our data.

We first find the most commonly expressed genes in our RNA data, then subset the vcf file to SNPs seen in those genes.





















