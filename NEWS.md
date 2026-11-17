# demuxSNP 1.4.0

Added reassign_centroid() function for improved handling of sparse and missing data.
Added messaging to previous functions reassign(), reassign_jaccard() and reassign_balanced() to point towards reassign_centroid().
Updated add_snps() function to add unfiltered data.

# demuxSNP 0.99

Submitted to Bioconductor

## Features

Functions to create a supervised demultiplexing classification model using 
  cell hashing and SNPs to aid in:  
* SNP selection
* creating training data using demuxmix
* adding SNPs data to SingleCellExperiment objects
* cell assignment


