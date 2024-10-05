
<!-- README.md is generated from README.Rmd. Please edit that file -->

# demuxSNP <img src="man/figures/logo.png" align="right" height="138" alt="" />

<!-- badges: start -->
<!-- badges: end -->

## Introduction

Multiplexing in scRNAseq involves the sequencing of samples from
different patients, treatment types or physiological locations together,
resulting in significant cost savings. The cells must then be
demultiplexed, or assigned back to their respective groups. A number of
experimental and computational methods have been proposed to facilitate
this, but a universally robust algorithm remains elusive. Below, we
introduce some existing methods, highlight the novel features of our
approach and its advantages to the user.

## Existing Methods

### Cell Hashing

Cells from each group are labelled with a distinct tag (HTO or LMO)
which is sequenced to give a counts matrix. Due to non-specific binding,
these counts form a bimodal distribution. Such methods are generally
computationally efficient. Their classification performance, however, is
highly dependent on the tagging quality and many methods do not account
for uncertainty in classification (Boggy et al. (2022), Kim et al.
(2020) & Stoeckius et al. (2018)).

More recent methods, including
[demuxmix](https://bioconductor.org/packages/release/bioc/html/demuxmix.html),
assign a probability that a cell is from a particular group, or made up
of multiple groups (doublet). This allows users to define a cut-off
threshold for the assignment confidence. Accounting for uncertainty is
an important feature for these types of algorithms. But, while they give
the user greater flexibility in determining which cells to keep, this
ultimately results in a trade-off between keeping cells which cannot be
confidently called or discarding them - due to issues with tagging
quality rather than RNA quality.

### SNPs

The second class of methods exploits natural genetic variation between
cells and so can only be used where the groups are genetically distinct.
Demuxlet (Kang et al. (2018)) uses genotype information from each group
to classify samples. This genotyping incurs additional experimental
cost. To address this, Souporcell (Heaton et al. (2020)) and Vireo
(Huang, McCarthy, and Stegle (2019)) among other methods were developed
to classify cells based on their SNPs in an unsupervised manner. Without
prior knowledge of the SNPs associated with each group, these
unsupervised methods may confuse groups with lower cell counts for other
signals in the data.

Demuxlet remains the standard often used to benchmark other methods but
its more widespread adoption has been limited by the requirement of
sample genotype information.

## demuxSNP Motivation

**With cell hashing, we can confidently demultiplex *some* but not *all*
cells. Using these high confidence cells, we can learn the SNPs
associated with each group. This SNP information can then be used to
assign remaining cells (which we could not confidently call using cell
hashing) to their most similar group based on their SNP profile.**

Novel features:

- Uses both cell hashing and SNP data. Current methods are limited to
  using one or the other.

- Selects SNPs based on being located in a gene expressed in a large
  proportion of cells to reduce noise, computational cost and increase
  interpretability.

Impact:

- Users can visually confirm validity (or lack thereof) of existing
  demultiplexing results in a tangible manner.
- Users can recover otherwise high quality cells which could not be
  confidently assigned using other methods.
- Cells from groups which are present in lower proportions may be
  classified better than with unsupervised SNP approaches.

Note: the approach used here differs from most SNP methods in that it is
supervised. We attain knowledge of which SNPs are associated with which
patients from the high confidence cells then use this to train a
classifier. It is similar to demuxlet in the sense that the classifier
uses group specific SNP information, **however** our method does not
require the expense of genotyping and so may be much more widely
applicable.

## Installation

You can install the development version of demuxSNP from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("michaelplynch/demuxSNP")
```

## Workflow

``` r
library(demuxSNP)
library(ComplexHeatmap)
library(viridisLite)
library(Seurat)
library(ggpubr)
library(dittoSeq)
library(utils)
library(EnsDb.Hsapiens.v86)
```

``` r
colors <- structure(viridis(n = 3), names = c("-1", "0", "1"))
```

``` r
#Data loading
data(multiplexed_scrnaseq_sce, 
     commonvariants_1kgenomes_subset, 
     vartrix_consensus_snps, 
     package = "demuxSNP")

small_sce<-multiplexed_scrnaseq_sce[,1:100]
ensdb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86

#Preprocessing
top_genes<-common_genes(small_sce)
vcf_sub<-subset_vcf(commonvariants_1kgenomes_subset, top_genes, ensdb)
small_sce<-high_conf_calls(small_sce)

#Use subsetted vcf with VarTrix in default 'consensus' mode to generate SNPs 
#matrix

small_sce<-add_snps(small_sce,vartrix_consensus_snps[,1:100])

small_sce<-reassign(small_sce,k=5)

table(small_sce$knn)
#> 
#>  Doublet Hashtag1 Hashtag3 Hashtag4 Hashtag5 Hashtag6 
#>       22       10       17        7       11       33
```

## Example

# Exploratory Analysis

We load three data objects: a SingleCellExperiment object containing RNA
and HTO counts, a .vcf file of class CollapsedVCF containing SNPs from
1000 Genomes common variants and a matrix containing SNP information for
each cell (we will show you how to generate this SNPs matrix using
[VarTrix](https://github.com/10XGenomics/vartrix) outside of R). We have
already removed low quality cells (library size\<1,000 and percentage of
genes mapping to mitochondrial genes\>10%).

``` r

class(multiplexed_scrnaseq_sce)
#> [1] "SingleCellExperiment"
#> attr(,"package")
#> [1] "SingleCellExperiment"
class(commonvariants_1kgenomes_subset)
#> [1] "CollapsedVCF"
#> attr(,"package")
#> [1] "VariantAnnotation"
class(vartrix_consensus_snps)
#> [1] "matrix" "array"
```

The HTO or LMO distribution is usually bimodal, with a signal (high
counts) and background distribution (low counts) caused by non-specific
binding. Ideally, these distributions would be clearly separated with no
overlap, but in practice, this is not always the case. In our example
data, we see that the signal and noise overlap to varying extents in
each group.

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" />

We will begin by running Seurat’s HTODemux, a popular HTO demultiplexing
algorithm on the data.

``` r
logcounts(multiplexed_scrnaseq_sce) <- counts(multiplexed_scrnaseq_sce)
seurat <- as.Seurat(multiplexed_scrnaseq_sce)
seurat <- HTODemux(seurat)
#> As of Seurat v5, we recommend using AggregateExpression to perform pseudo-bulk analysis.
#> First group.by variable `ident` starts with a number, appending `g` to ensure valid variable names
#> Cutoff for Hashtag1 : 44 reads
#> 
#> Cutoff for Hashtag2 : 39 reads
#> 
#> Cutoff for Hashtag3 : 323 reads
#> 
#> Cutoff for Hashtag4 : 99 reads
#> 
#> Cutoff for Hashtag5 : 107 reads
#> 
#> Cutoff for Hashtag6 : 175 reads
#> 
#> This message is displayed once per session.
seurat$hash.ID <- factor(as.character(seurat$hash.ID))
multiplexed_scrnaseq_sce$seurat <- seurat$hash.ID

multiplexed_scrnaseq_sce$seurat <- seurat$hash.ID

table(multiplexed_scrnaseq_sce$seurat)
#> 
#>  Doublet Hashtag1 Hashtag2 Hashtag3 Hashtag4 Hashtag5 Hashtag6 Negative 
#>      633      121       29      264      158      177      383      235
```

Although HTO library size of the Negative group is low, the RNA library
size is similar to that of other groups, indicating that they may be
misclassified as Negative due to their tagging quality rather than
overall RNA quality.

``` r
dittoPlot(seurat, "nCount_HTO", group.by = "hash.ID")
```

<img src="man/figures/README-unnamed-chunk-9-1.png" width="100%" />

``` r
dittoPlot(seurat, "nCount_RNA", group.by = "hash.ID")
```

<img src="man/figures/README-unnamed-chunk-9-2.png" width="100%" />

For the remainder of this vignette we will outline our method of
checking whether or not these cells have been called correctly and how
to assign them to their appropriate group!

# Preprocessing

Common variants files, for example from the 1000 Genomes Project, can
contain over 7 million SNPs. To reduce computational cost and cell-type
effects, we subset our SNPs list to those located within genes expressed
across most cells in our data.

We first find the most commonly expressed genes in our RNA data.

``` r
top_genes <- common_genes(sce = multiplexed_scrnaseq_sce, n = 100)

top_genes[1:10]
#>  [1] "TPT1"   "RPL13"  "RPL28"  "TMSB4X" "RPS27"  "EEF1A1" "RPL41"  "B2M"   
#>  [9] "RPLP1"  "RPL32"
```

We have a sample .vcf preloaded, but you can load your .vcf file in
using ‘readVcf()’ from
[VariantAnnotation](https://bioconductor.org/packages/release/bioc/html/VariantAnnotation.html).

We will subset our .vcf file to SNPs seen in commonly expressed genes
from our dataset. Notice that the genome for the vcf and EnsDb object
must be compatible!

The returned vcf can be written to file and used with VarTrix later on.

``` r
ensdb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86

genome(commonvariants_1kgenomes_subset)[1] == genome(ensdb)[1]
#>    1 
#> TRUE

new_vcf <- subset_vcf(commonvariants_1kgenomes_subset, top_genes = top_genes, ensdb)
commonvariants_1kgenomes_subset
#> class: CollapsedVCF 
#> dim: 2609 0 
#> rowRanges(vcf):
#>   GRanges with 5 metadata columns: paramRangeID, REF, ALT, QUAL, FILTER
#> info(vcf):
#>   DataFrame with 1 column: AF
#> info(header(vcf)):
#>       Number Type  Description                                  
#>    AF A      Float Estimated allele frequency in the range (0,1)
#> geno(vcf):
#>   List of length 0:
new_vcf
#> class: CollapsedVCF 
#> dim: 2399 0 
#> rowRanges(vcf):
#>   GRanges with 5 metadata columns: paramRangeID, REF, ALT, QUAL, FILTER
#> info(vcf):
#>   DataFrame with 1 column: AF
#> info(header(vcf)):
#>       Number Type  Description                                  
#>    AF A      Float Estimated allele frequency in the range (0,1)
#> geno(vcf):
#>   List of length 0:
```

The subsetted .vcf can be written to disk using ‘writeVcf()’, again from
[VariantAnnotation](https://bioconductor.org/packages/release/bioc/html/VariantAnnotation.html)
package.

Next, we wish to identify cells which we can confidently call to a
particular group. There are a number of ways this can be achieved,
including probabilistic modelling of the HTO counts, manually setting
non-conservative thresholds or using consensus calls. The user may wish
to experiment with different approaches. Here we have used
[demuxmix](https://bioconductor.org/packages/release/bioc/html/demuxmix.html),
a probabilistic model which we have set with a high acceptance threshold
to identify high confidence cell calls to use as training data (cells
which we can confidently call as a particular singlet group).

``` r
multiplexed_scrnaseq_sce <- high_conf_calls(multiplexed_scrnaseq_sce)

table(multiplexed_scrnaseq_sce$train)
#> 
#> FALSE  TRUE 
#>   955  1045

table(multiplexed_scrnaseq_sce$predict)
#> 
#> TRUE 
#> 2000

table(multiplexed_scrnaseq_sce$labels)
#> 
#>  Hashtag1  Hashtag2  Hashtag3  Hashtag4  Hashtag5  Hashtag6 multiplet  negative 
#>        62        15       226       102       348       292       335        12 
#> uncertain 
#>       608
```

So, for this particular dataset, we can confidently call 1,045 cells as
being from a particular singlet group. 608 cells cannot be called to a
group with high confidence.

# Variant Calling (VarTrix)

Variant calling is not done within the package. Instead, we refer the
reader to [VarTrix](https://github.com/10XGenomics/vartrix), where they
can use the subsetted .vcf file along with their .bam, barcodes.tsv and
reference genome to call SNPs in each cell.

A sample VarTrix command looks like the following:

``` bash

./vartrix -v <path_to_input_vcf> -b <path_to_cellranger_bam> -f <path_to_fasta_file> -c <path_to_cell_barcodes_file> -o <path_for_output_matrix>
```

Using the output matrix from Vartrix and the high confidence
classifications from the HTO algorithm, we can reassign cells using
k-nearest neighbours.

# Cell Reassignment, Visualisation and Evaluation

To keep things tidy, we will add the SNP data to our
SingleCellExperiment object as an altExp. We recode the SNP matrix as
follows: 0=no read, 1=SNP present, -1=SNP absent. This function also
filters out SNPs which are observed at a low frequency in the data, and
the frequency threshold can be set manually.

``` r
dim(vartrix_consensus_snps)
#> [1] 2542 2000

multiplexed_scrnaseq_sce <- add_snps(multiplexed_scrnaseq_sce, vartrix_consensus_snps, thresh = 0.95)

altExp(multiplexed_scrnaseq_sce, "SNP")
#> class: SingleCellExperiment 
#> dim: 85 2000 
#> metadata(0):
#> assays(1): counts
#> rownames(85): Snp Snp ... Snp Snp
#> rowData names(0):
#> colnames(2000): AAACCTGAGATCTGCT-1 AAACCTGAGCGTCAAG-1 ...
#>   ACTTTCAGTAAGTTCC-1 ACTTTCAGTAGGCATG-1
#> colData names(0):
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
```

Before we reassign any cells, we will first use the SNPs data to inspect
the results from stand-alone algorithms. Splitting the SNP data by
Seurat HTODemux classification, we initially see a large number of
‘negative’ cells which appear of good quality (high proportion of reads)
which may be assignable to another group. This is consistent with the
library size plot we visualised earlier.

``` r
hm <- Heatmap(counts(altExp(multiplexed_scrnaseq_sce, "SNP")),
    column_split = multiplexed_scrnaseq_sce$seurat,
    cluster_rows = FALSE,
    show_column_names = FALSE,
    cluster_column_slices = FALSE,
    column_title_rot = -45,
    row_title = "SNPs",
    show_row_names = FALSE,
    col = colors
)

draw(hm,
    column_title = "SNP profiles split by Seurat HTODemux call",
    padding = unit(c(2, 15, 2, 2), "mm")
)
```

<img src="man/figures/README-unnamed-chunk-15-1.png" width="100%" />

We will use our knn method to reassign cells based on their SNP
profiles. The training data is the high confidence cells

``` r
set.seed(1)
multiplexed_scrnaseq_sce$labels<-as.character(multiplexed_scrnaseq_sce$labels)
multiplexed_scrnaseq_sce$knn <- reassign_centroid(multiplexed_scrnaseq_sce,
    train_cells = multiplexed_scrnaseq_sce$train,
    predict_cells = multiplexed_scrnaseq_sce$predict,
    min_cells = 5
)
#> [1] 2542    6
#> [1] 3837
#> Hashtag1 Hashtag2 Hashtag3 Hashtag4 Hashtag5 Hashtag6  Doublet  Doublet 
#>      579      394      845      571      696      752      501      756 
#>  Doublet  Doublet  Doublet  Doublet  Doublet  Doublet  Doublet  Doublet 
#>      689      746      743      528      480      499      477      773 
#>  Doublet  Doublet  Doublet  Doublet  Doublet 
#>      875      944      706      719      846 
#>                    Hashtag1 Hashtag2 Hashtag3 Hashtag4 Hashtag5 Hashtag6
#> AAACCTGAGATCTGCT-1       71       47       60       53       52       51
#> AAACCTGAGCGTCAAG-1       46       48       40       61       44       49
#> AAACCTGAGGCGTACA-1       40       49       49       45       51       68
#> AAACCTGAGGGCTCTC-1       64       78       61       59       62       73
#> AAACCTGAGTAGGTGC-1       63       65       65       64       73       91
#>                    Hashtag1 Hashtag2 Hashtag3 Hashtag4 Hashtag5 Hashtag6
#> AAACCTGAGATCTGCT-1        7       23       19       25       26       27
#> AAACCTGAGCGTCAAG-1       16       12       24        3       20       15
#> AAACCTGAGGCGTACA-1       27       14       19       22       15        0
#> AAACCTGAGGGCTCTC-1       36        6       45       42       41       31
#> AAACCTGAGTAGGTGC-1       36       19       34       33       25        8
#>                    Hashtag1 Hashtag2 Hashtag3 Hashtag4 Hashtag5 Hashtag6
#> AAACCTGAGATCTGCT-1        3       23       20       17       24       22
#> AAACCTGAGCGTCAAG-1       20       19       36        8       26       23
#> AAACCTGAGGCGTACA-1       29       19       26       23       17        3
#> AAACCTGAGGGCTCTC-1       23        4       36       23       29       25
#> AAACCTGAGTAGGTGC-1       31       21       29       25       20        4
#>                    Hashtag1 Hashtag2 Hashtag3 Hashtag4 Hashtag5 Hashtag6
#> AAACCTGAGATCTGCT-1      109       75       94       95       90       91
#> AAACCTGAGCGTCAAG-1       96       84       82      109       92       94
#> AAACCTGAGGCGTACA-1       85       82       94       91       99      118
#> AAACCTGAGGGCTCTC-1      118      108      111      119      119      120
#> AAACCTGAGTAGGTGC-1      113      102      115      116      122      143
#> [1] 1500

table(multiplexed_scrnaseq_sce$knn)
#> 
#>  Doublet Hashtag1 Hashtag2 Hashtag3 Hashtag4 Hashtag5 Hashtag6 
#>      464      119       30      340      155      429      463
```

``` r
hm <- Heatmap(counts(altExp(multiplexed_scrnaseq_sce, "SNP")),
    column_split = multiplexed_scrnaseq_sce$knn,
    cluster_rows = FALSE,
    show_column_names = FALSE,
    cluster_column_slices = FALSE,
    column_names_rot = 45,
    column_title_rot = -45,
    row_title = "SNPs",
    show_row_names = FALSE,
    col = colors
)

draw(hm,
    column_title = "SNP profiles split by updated knn classification",
    padding = unit(c(2, 15, 2, 2), "mm")
)
```

<img src="man/figures/README-unnamed-chunk-17-1.png" width="100%" />

Focusing in on the knn Hashtag5 group, we see that a lot of the Negative
cells have now been correctly reclassed to this group, as well as a
small number of cells from other groups.

``` r
hm <- Heatmap(counts(altExp(multiplexed_scrnaseq_sce, "SNP"))[, multiplexed_scrnaseq_sce$knn == "Hashtag5"],
    column_split = multiplexed_scrnaseq_sce$seurat[multiplexed_scrnaseq_sce$knn == "Hashtag5"],
    cluster_rows = FALSE,
    show_column_names = FALSE,
    cluster_column_slices = FALSE,
    column_names_rot = 45,
    column_title_rot = -45,
    row_title = "SNPs",
    show_row_names = FALSE,
    col = colors
)

draw(hm,
    column_title = "knn Hashtag5 group split by Seurat HTODemux classification",
    padding = unit(c(2, 15, 2, 2), "mm")
)
```

<img src="man/figures/README-unnamed-chunk-18-1.png" width="100%" />

## Performance

Next we will run some basic performance checks. We subset our
SingleCellExperiment object to only retain cells which we could
confidently call as being from a singlet group, then split this into a
training and test dataset.

``` r
sce_test <- multiplexed_scrnaseq_sce[, multiplexed_scrnaseq_sce$train == TRUE]
sce_test$knn <- NULL
#sce_test$labels <- droplevels(sce_test$labels)
sce_test
#> class: SingleCellExperiment 
#> dim: 259 1045 
#> metadata(0):
#> assays(2): counts logcounts
#> rownames(259): RPL22 CDC42 ... MT-ND5 MT-CYB
#> rowData names(0):
#> colnames(1045): AAACCTGAGCGTCAAG-1 AAACCTGAGGCGTACA-1 ...
#>   ACTTTCAGTAAGAGAG-1 ACTTTCAGTAAGTTCC-1
#> colData names(11): orig.ident nCount_RNA ... predict labels
#> reducedDimNames(0):
#> mainExpName: RNA
#> altExpNames(3): HTO SNPcons SNP

sce_test$train2 <- rep(FALSE, length(sce_test$train))
sce_test$train2[seq_len(500)] <- TRUE

sce_test$test <- sce_test$train2 == FALSE
```

Comparing the predicted labels in the test dataset with the hidden high
confidence labels, we see excellent agreement.

``` r
sce_test$knn <- reassign_centroid(sce_test, train_cells = sce_test$train2, predict_cells = sce_test$test, min_cells = 5)
#> [1] 2542    6
#> [1] 3837
#> Hashtag1 Hashtag2 Hashtag3 Hashtag4 Hashtag5 Hashtag6  Doublet  Doublet 
#>      579      394      845      571      696      752      501      756 
#>  Doublet  Doublet  Doublet  Doublet  Doublet  Doublet  Doublet  Doublet 
#>      689      746      743      528      480      499      477      773 
#>  Doublet  Doublet  Doublet  Doublet  Doublet 
#>      875      944      706      719      846 
#>                    Hashtag1 Hashtag2 Hashtag3 Hashtag4 Hashtag5 Hashtag6
#> AAACCTGAGCGTCAAG-1       46       48       40       61       44       49
#> AAACCTGAGGCGTACA-1       40       49       49       45       51       68
#> AAACCTGAGGGCTCTC-1       64       78       61       59       62       73
#> AAACCTGAGTAGGTGC-1       63       65       65       64       73       91
#> AAACCTGAGTCAATAG-1       56       58       58       61       66       83
#>                    Hashtag1 Hashtag2 Hashtag3 Hashtag4 Hashtag5 Hashtag6
#> AAACCTGAGCGTCAAG-1       16       12       24        3       20       15
#> AAACCTGAGGCGTACA-1       27       14       19       22       15        0
#> AAACCTGAGGGCTCTC-1       36        6       45       42       41       31
#> AAACCTGAGTAGGTGC-1       36       19       34       33       25        8
#> AAACCTGAGTCAATAG-1       28       16       33       23       19        7
#>                    Hashtag1 Hashtag2 Hashtag3 Hashtag4 Hashtag5 Hashtag6
#> AAACCTGAGCGTCAAG-1       20       19       36        8       26       23
#> AAACCTGAGGCGTACA-1       29       19       26       23       17        3
#> AAACCTGAGGGCTCTC-1       23        4       36       23       29       25
#> AAACCTGAGTAGGTGC-1       31       21       29       25       20        4
#> AAACCTGAGTCAATAG-1       28       19       30       27       23        4
#>                    Hashtag1 Hashtag2 Hashtag3 Hashtag4 Hashtag5 Hashtag6
#> AAACCTGAGCGTCAAG-1       96       84       82      109       92       94
#> AAACCTGAGGCGTACA-1       85       82       94       91       99      118
#> AAACCTGAGGGCTCTC-1      118      108      111      119      119      120
#> AAACCTGAGTAGGTGC-1      113      102      115      116      122      143
#> AAACCTGAGTCAATAG-1      122      110      130      126      136      159

table(sce_test$labels[sce_test$test == TRUE], sce_test$knn[sce_test$test == TRUE])
#>           
#>            Doublet Hashtag1 Hashtag2 Hashtag3 Hashtag4 Hashtag5 Hashtag6
#>   Hashtag1       2       34        0        0        0        0        0
#>   Hashtag2       0        0       11        0        0        0        0
#>   Hashtag3       0        0        0      122        0        0        0
#>   Hashtag4       1        0        0        0       45        0        0
#>   Hashtag5       7        0        0        0        0      187        1
#>   Hashtag6       1        0        0        0        0        0      134
```

# References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-boggy_bff_2022" class="csl-entry">

Boggy, Gregory J, G W McElfresh, Eisa Mahyari, Abigail B Ventura, Scott
G Hansen, Louis J Picker, and Benjamin N Bimber. 2022. “BFF and
<span class="nocase">cellhashR</span>: Analysis Tools for Accurate
Demultiplexing of Cell Hashing Data.” *Bioinformatics* 38 (10):
2791–801. <https://doi.org/10.1093/bioinformatics/btac213>.

</div>

<div id="ref-heaton_souporcell_2020" class="csl-entry">

Heaton, Haynes, Arthur M. Talman, Andrew Knights, Maria Imaz, Daniel J.
Gaffney, Richard Durbin, Martin Hemberg, and Mara K. N. Lawniczak. 2020.
“Souporcell: Robust Clustering of Single-Cell RNA-Seq Data by Genotype
Without Reference Genotypes.” *Nature Methods* 17 (6): 615–20.
<https://doi.org/10.1038/s41592-020-0820-1>.

</div>

<div id="ref-huang_vireo_2019" class="csl-entry">

Huang, Yuanhua, Davis J. McCarthy, and Oliver Stegle. 2019. “Vireo:
Bayesian Demultiplexing of Pooled Single-Cell RNA-Seq Data Without
Genotype Reference.” *Genome Biology* 20 (1): 273.
<https://doi.org/10.1186/s13059-019-1865-2>.

</div>

<div id="ref-kang_multiplexed_2018" class="csl-entry">

Kang, Hyun Min, Meena Subramaniam, Sasha Targ, Michelle Nguyen, Lenka
Maliskova, Elizabeth McCarthy, Eunice Wan, et al. 2018. “Multiplexed
Droplet Single-Cell RNA-Sequencing Using Natural Genetic Variation.”
*Nature Biotechnology* 36 (1): 89–94.
<https://doi.org/10.1038/nbt.4042>.

</div>

<div id="ref-kim_citefuse_2020" class="csl-entry">

Kim, Hani Jieun, Yingxin Lin, Thomas A Geddes, Jean Yee Hwa Yang, and
Pengyi Yang. 2020. “CiteFuse Enables Multi-Modal Analysis of CITE-Seq
Data.” *Bioinformatics* 36 (14): 4137–43.
<https://doi.org/10.1093/bioinformatics/btaa282>.

</div>

<div id="ref-stoeckius_cell_2018" class="csl-entry">

Stoeckius, Marlon, Shiwei Zheng, Brian Houck-Loomis, Stephanie Hao,
Bertrand Z. Yeung, William M. Mauck, Peter Smibert, and Rahul Satija.
2018. “Cell Hashing with Barcoded Antibodies Enables Multiplexing and
Doublet Detection for Single Cell Genomics.” *Genome Biology* 19 (1):
224. <https://doi.org/10.1186/s13059-018-1603-1>.

</div>

</div>
