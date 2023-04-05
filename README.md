# demuxSNP

R package for demultiplexing scRNAseq data. Uses cell hashing and SNPs to QC and reassign miscalled or uncalled cells.

Multiplexing in scRNAseq involves pooling and loading multiple biological samples onto the same sequencing lane to reduce sequencing costs.
Demultiplexing must then be  carried out to identify which cells came from which biological sample.

Demultiplexing methods fall broadly into two categories, either using cell hashing or SNPs.
demuxSNP is of particular use in applications where cell hashing quality is poor (impacting performance of hashing based methods) and some sample groups may have low cell numbers (impacting performance of genotype free SNPs methods).

# Install

To install the development version:

```
devtools::install_github("michaelplynch/demuxSNP")
```

# How to

The demuxSNP package contains functions to assist the user in key parts of the workflow.

1. `common_genes` and `subset_vcf` find widely observed genes in the dataset 
then subset a SNPs vcf file to SNPs seen within those genes.
2. `high_conf_calls` wraps around demuxmix and marks singlet cells assigned 
with high confidence to use to train the model.
3. `add_snps` adds the snps matrix (default 'consensus' output from [VarTrix](https://github.com/10XGenomics/vartrix)) to the SingleCellExperiment 
object and carries out additional filtering to remove SNPs not observed in the 
data.
4. `reassign` trains a knn classifier using singlets in the training data 
(from high_conf_calls function) and simulated doublets, and uses this to 
classify cells based on their SNP profile.


Please see package vignette for a full demonstration of functions.
