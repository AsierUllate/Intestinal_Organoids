####
# Script to obtain R object for human Intestine crypt from Fujii paper
###
library(Seurat)
library(tidyverse)
library(findPC)
library(future)
plan("multiprocess", workers = 6)

#Data in GSM3389578
#LogNormalized counts for intestine, have only cells after filtering
logcounts_Intestine <- read.delim("data/Fujii_2018/GSM3389578_tissue_LogNormalized.txt",
sep = "\t")

#Raw counts for intestine, have all cells
counts_Intestine <- read.delim("data/Fujii_2018/GSM3389578_tissue_count.txt",
sep = "\t")

rownames(counts_Intestine) <- counts_Intestine$Gene
counts_Intestine <- counts_Intestine[,-1]

#Keep only cells that are there after filtering
counts_Intestine <- counts_Intestine[, colnames(counts_Intestine) %in% colnames((logcounts_Intestine))]

colnames(counts_Intestine) <- str_replace(colnames(counts_Intestine), pattern = "\\.", replacement = "-")

#Check if RANKL, RANK and OPG genes are in the counts matrix
genes <- c("TNFRSF11A", "TNFSF11")
intersect(genes, rownames(counts_Intestine)) #RANKL is there as a gene

Intestine <- CreateSeuratObject(counts = counts_Intestine, project = "Intestine_Organoid")
Intestine$percent.mt <- PercentageFeatureSet(Intestine, pattern = "^MT-")

Intestine <- NormalizeData(Intestine) %>%
FindVariableFeatures() %>%
ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA")) %>%
RunPCA()

Intestine <- RunUMAP(Intestine, dims = 1:10)

Intestine <- FindNeighbors(Intestine, dims = 1:10)
Intestine <- FindClusters(Intestine, res = 0.6)

p1 <- DimPlot(Intestine)
ggratio <- diff(range(p1$data[, 2])) / diff(range(p1$data[, 1]))

#Check how it looks at different resolutions
cols_clusters <- RColorBrewer::brewer.pal("Paired", n = 12)

pdf("UMAPs_resolutions_Intestine.pdf")
p1 <- DimPlot(Intestine, label = T, cols = c(cols_clusters), group.by = "RNA_snn_res.0.6") & coord_fixed(ratio = ggratio)
print(p1)
dev.off()

#They used 25 PCs in the paper, with k.param=20 in Intestine for Findneighbors and resolution 0.5
# Perplexity of 30 in RuntsNE
Intestine_Paper <- FindNeighbors(Intestine, dims = 1:25, k.param = 20) %>%
FindClusters(res = 0.6)

Intestine_Paper <- RunUMAP(Intestine_Paper, dims = 1:25)

p1 <- DimPlot(Intestine_Paper)
ggratio <- diff(range(p1$data[, 2])) / diff(range(p1$data[, 1]))

pdf(file.path(project_dir, "UMAPs_resolutions_Intestine_25PCs.pdf"))

DimPlot(Intestine_Paper, label = T, cols = cols_clusters, group.by = "RNA_snn_res.0.6") & coord_fixed(ratio = ggratio)

dev.off()

#Remove cluster 9 at resolution 0.6 from Intestine as non-epithelial
Intestine <- subset(Intestine, subset = RNA_snn_res.0.6 != "9")

Intestine <- FindVariableFeatures(Intestine) %>%
ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA")) %>%
RunPCA()

#They used 25 PCs in the paper, with k.param=40 in Intestine for Findneighbors and resolution 0.5
# Perplexity of 30 in RuntsNE
Intestine_Paper <- FindNeighbors(Intestine, dims = 1:25, k.param = 20) %>%
FindClusters(res = 0.8)

Intestine_Paper <- RunUMAP(Intestine_Paper, dims = 1:25)

p1 <- DimPlot(Intestine_Paper)
ggratio <- diff(range(p1$data[, 2])) / diff(range(p1$data[, 1]))

pdf("UMAPs_resolutions_Intestine_25PCs_epithelialonly.pdf")
p1 <- DimPlot(Intestine_Paper, label = T, cols = c(cols_clusters, "lightgray"), group.by = "RNA_snn_res.0.8") & coord_fixed(ratio = ggratio)
print(p1)
dev.off()


Idents(Intestine_Paper) <- "RNA_snn_res.0.8"
DefaultAssay(Intestine_Paper) <- "RNA"
Intestine_Paper <- NormalizeData(Intestine_Paper)
#Get markers for each of the clusters
# annotations <- read_tsv("/home/aullatea/Penninger/annotation_Mus_v98.tsv")
temp<- FindAllMarkers(Intestine_Paper, min.pct = 0, logfc.threshold = 0, only.pos = T) %>% mutate(specific = pct.1 - pct.2)
temp <- temp %>% relocate(cluster) %>% arrange(cluster, desc(avg_log2FC)) %>% filter(avg_log2FC > 0.25 & specific > 0.1)

write_tsv(temp, file = file.path(project_dir, "Markers_ClustersIntestine_res0.8_DEG.xls"))

#Put those that we consider with the same annotation together
Intestine_Paper$CellType <- factor(Intestine_Paper$RNA_snn_res.0.8,
levels = levels(Intestine_Paper$RNA_snn_res.0.8),
labels = c("Stem cells", "Stem cells", "TA 1 cells", "Early Enterocytes", "Goblet", "TA 1 cells", "TA 2 cells", "TA 2 cells",
"Enterocytes", "Paneth cells", "Stem cells", "Microfold cells", "Tuft/Early EECs"))

Idents(Intestine_Paper) <- "CellType"
Intestine_Paper$RNA_snn_res.0.6 <- NULL
saveRDS(Intestine_Paper, "data/Fujii_2018/Fujii.rds")

