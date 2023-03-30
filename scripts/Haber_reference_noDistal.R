####
# Script used to generate the Haber reference
####
library(Seurat)
library(tidyverse)
library(ggpubr)
library(future)
options(future.globals.maxSize = 12000 * 1024^2)
options(bitmapType = "cairo")
plan("multiprocess", workers = 8)

# Already normalized count matrix from droplet-based scRNAseq from 6 mice and metadata
# Obtained from https://singlecell.broadinstitute.org/single_cell/study/SCP44/small-intestinal-epithelium#study-download
counts <- read.table("atlas_Log2Tpm_round2.txt", row.names = 1, header = T, sep = "\t", as.is = T)
meta <- read_tsv("atlas_metadata.txt")

# Get the UMI counts from NCBI GEO GSE92332
counts2 <- read.table("GSE92332_atlas_UMIcounts.txt", row.names = 1, header = T, sep = "\t", as.is = T)

# Compare barcodes to check that the order is the same
barcodes <- sapply(str_split(colnames(counts), pattern = "_"), "[[", 1)
barcodes2 <- sapply(str_split(colnames(counts2), pattern = "_"), "[[", 2)
all(barcodes == barcodes2)
samp <- sapply(str_split(colnames(counts2), pattern = "_"), "[[", 1)
celltype <- sapply(str_split(colnames(counts2), pattern = "_"), "[[", 3)
# Change the barcodes name from the counts2 table to the one from the log2 to use the metadata
colnames(counts2) <- colnames(counts)


meta <- meta[-1, ]
meta <- meta[meta$NAME %in% colnames(counts), ]
meta <- as.data.frame(meta)
rownames(meta) <- meta$NAME
meta <- meta[colnames(counts), ]
meta$batch <- samp

# Check how different are the gene names from the database compared to the one used
genes_2020 <- read_tsv("features_2020.tsv", col_names = F)
genes_old <- read_tsv("GSM2839445_Atlas1_genes.tsv", col_names = F)
gene_names <- Read10X_h5("/home/aullatea/Penninger/scRNA_WT_Organoid/scRNA_WT_Organoid_Count/outs/filtered_feature_bc_matrix.h5")

# Get actual gene names as for our analyses
genes_2020$X2 <- rownames(gene_names)
diff <- setdiff(genes_old$X1, genes_2020$X1) # 979 Ensembl Ids are not in the latest version anymore
nonexisting <- intersect(genes_old$X2[which(genes_old$X1 %in% diff)], rownames(counts))
length(nonexisting) # 352 genes, 2.2% of genes
sum(rowSums(counts2[nonexisting, ])) # 108490
sum(rowSums(counts2[nonexisting, ])) * 100 / sum(counts2) # 0.2% of UMIs
# Convert to the same symbols as for us
genes_old$new <- sapply(genes_old$X1, function(x) {
  ifelse(x %in% genes_2020$X1, genes_2020$X2[which(genes_2020$X1 == x)], genes_old$X2[which(genes_old$X1 == x)])
})

# Change rownames to Ensembl IDs to avoid issues with double names
rownames(counts2) <- sapply(rownames(counts2), function(x) {
  ifelse(x %in% genes_old$X2, genes_old$X1[which(genes_old$X2 == x)], str_split(x, pattern = "_")[[1]][2])
})
# Remove those that are not in the latest version
counts2 <- counts2[!(rownames(counts2) %in% diff), ]
rownames(counts2) <- genes_old$new[match(rownames(counts2), genes_old$X1)]

# Create a Seurat object
Sample1 <- CreateSeuratObject(counts = counts2, project = "Haber", meta.data = meta, min.cells = 3, min.features = 200)

# Remove distal cells
Sample1 <- subset(Sample1, subset = Cluster != "Enterocyte.Mature.Distal" & Cluster != "Enterocyte.Immature.Distal")
Sample1$Cluster_SimplifiedTA <- factor(Sample1$Cluster,
  levels = c(
    "Enterocyte.Immature.Proximal", "Enterocyte.Mature.Proximal",
    "Enterocyte.Progenitor", "Enterocyte.Progenitor.Early",
    "Enterocyte.Progenitor.Late", "Enteroendocrine", "Goblet",
    "Paneth", "Stem", "TA", "TA.G1", "TA.G2", "Tuft"
  ),
  labels = c(
    "Enterocyte.Immature.Proximal", "Enterocyte.Mature.Proximal",
    "Enterocyte.Progenitor", "Enterocyte.Progenitor.Early",
    "Enterocyte.Progenitor.Late", "Enteroendocrine", "Goblet",
    "Paneth", "Stem", "TA", "TA", "TA", "Tuft"
  )
)
# it is mouse, so the names are in lowercase. There are no mitochondrial genes
Sample1[["percent.mt"]] <- PercentageFeatureSet(Sample1, pattern = "^mt-")
Sample1$percent.ribo <- PercentageFeatureSet(Sample1, pattern = "^Rp[sl]")

Sample1 <- NormalizeData(Sample1)
load("mouse_CC.Rdata")

Sample1 <- CellCycleScoring(Sample1, s.features = m.s.genes, g2m.features = m.g2m.genes)
Sample1$CC.Difference <- Sample1$S.Score - Sample1$G2M.Score


list_samples <- SplitObject(Sample1, split.by = "Mouse")
list_samples <- lapply(list_samples, SCTransform)

# Select integration features
samples.features <- SelectIntegrationFeatures(object.list = list_samples, nfeatures = 3000)
list_samples <- PrepSCTIntegration(
  object.list = list_samples, anchor.features = samples.features,
  verbose = FALSE
)

samples.anchors <- FindIntegrationAnchors(object.list = list_samples, normalization.method = "SCT", anchor.features = samples.features, verbose = TRUE)
samples.integrated <- IntegrateData(anchorset = samples.anchors, normalization.method = "SCT", verbose = TRUE)

samples.integrated <- RunPCA(samples.integrated, verbose = FALSE)
samples.integrated <- RunUMAP(samples.integrated, dims = 1:20)
samples.integrated <- RunTSNE(samples.integrated, dims = 1:20)

# Save the reference to label our samples later on
saveRDS(samples.integrated, "Haber_integrated_noDistal.rds")