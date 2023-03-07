## code to prepare `DATASET` dataset goes here
library(Matrix)
library(Seurat)
library(SingleCellExperiment)
## sce object

file_path <- c("C:/Users/michael.lynch/Culhane_Lab Dropbox/Shared_Lab_Folder/broad_datasets/KW9275_Yufei/210827_10X_KW9275_bcl/cellranger-6.1.1/GRCh38/BRI-1348/outs/filtered_feature_bc_matrix")
dat <- Read10X(file_path)

seurat <- CreateSeuratObject(dat$`Gene Expression`)
DefaultAssay(seurat) <- "RNA"
seurat <- NormalizeData(seurat) %>%
    FindVariableFeatures() %>%
    ScaleData()

seurat[["HTO"]] <- CreateAssayObject(counts = dat$`Antibody Capture`[1:6, ])
seurat <- NormalizeData(seurat, assay = "HTO", normalization.method = "CLR") %>% ScaleData()

seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")

seurat$qual <- rep("pass", each = 20700)
seurat$qual[seurat$percent.mt > 10] <- "fail"
seurat$qual[seurat$nCount_RNA < 1000] <- "fail"

sce <- as.SingleCellExperiment(seurat)
sce$ident <- NULL

keep <- sce$qual == "pass"

sce <- sce[, keep]

sce_sub <- sce[rowSums(counts(sce) > 0) > 14000, ]

logcounts(sce_sub) <- NULL
sce_sub <- sce_sub[, 1:2000]
sce <- sce_sub
usethis::use_data(sce, overwrite = TRUE)

# snps
filepath2 <- "C:/Users/michael.lynch/Culhane_Lab Dropbox/Shared_Lab_Folder/SNPs/v_outs_full.mtx"
snps_full <- as.matrix(readMM(filepath2))

snps_full <- snps_full[, keep]
snps <- snps_full[, 1:2000]

usethis::use_data(snps, overwrite = TRUE)

# vcf
filepath3 <- "C:/Users/michael.lynch/Culhane_Lab Dropbox/Shared_Lab_Folder/SNPs/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.chr.vcf"
my_vcf <- readVcf(filepath3, genome = "GRCh38")
SNP_ranges <- SummarizedExperiment::rowRanges(my_vcf)

my_vcf_inbound <- my_vcf[BiocGenerics::end(SNP_ranges) <= GenomeInfoDb::seqlengths(SNP_ranges)[as.character(GenomeInfoDb::seqnames(SNP_ranges))]]
SNP_ranges_inbound <- SummarizedExperiment::rowRanges(my_vcf_inbound)

gns <- ensembldb::genes(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)

top_gene_ranges <- gns[gns$gene_name %in% top_genes]

GenomeInfoDb::seqlengths(SNP_ranges_inbound) <- NA
GenomeInfoDb::seqlengths(top_gene_ranges) <- NA

include <- IRanges::overlapsAny(SNP_ranges_inbound, top_gene_ranges, type = "within")
include[sample(seq_len(14000), 250)] <- TRUE
top_genes_vcf <- my_vcf_inbound[include]
vcf <- top_genes_vcf

usethis::use_data(vcf, overwrite = FALSE)
