####
# Script to generate and annotate the integrated object
####

library(Seurat)
library(tidyverse)
library(future)
options(future.globals.maxSize = 12000 * 1024^2)
library(BiocParallel)
options(bitmapType = "cairo")
plan("multiprocess", workers = 8)


# Load CellCycle data
load("data/mouse_CC.Rdata")

wt <- Read10X_h5("/home/aullatea/Penninger/scRNA_WT_Organoid/scRNA_WT_Organoid_Count/outs/filtered_feature_bc_matrix.h5")
rankl <- Read10X_h5("/home/aullatea/Penninger/scRNA_RANKL_Organoid/scRNA_RANKL_Organoid_Count/outs/filtered_feature_bc_matrix.h5")

# Generate Seurat object for RANKL
rankl <- CreateSeuratObject(counts = rankl, project = "Organoids", min.cells = 3, min.features = 0)
rankl$Condition <- "RANKL"
rankl$orig.ident <- "RANKL_Organoid"

# Generate Seurat object for WT
wt <- CreateSeuratObject(counts = wt, project = "Organoids", min.cells = 3, min.features = 0)
wt$Condition <- "WT"
wt$orig.ident <- "WT_Organoid"

list_samples <- list("RANKL" = rankl, "Ctrl" = wt)

list_samples <- lapply(list_samples, function(x) {
  # it is mouse, so the names are in lowercase
  x$percent.mt <- PercentageFeatureSet(x, pattern = "^mt-")
  x$percent.ribo <- PercentageFeatureSet(x, pattern = "^Rp[sl]")
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
  x <- CellCycleScoring(x, s.features = m.s.genes, g2m.features = m.g2m.genes)
  x$CC.Difference <- x$S.Score - x$G2M.Score
  # Subset considering a minimun of 1,000 features and maximum of 20% MT
  x <- subset(x, subset = percent.mt < 20 & nFeature_RNA > 1000)
  x <- RenameCells(x, add.cell.id = x$orig.ident)
})

# Merge the two samples, to remove genes with low count levels and check how is the sample effect
samples.merged <- merge(x = list_samples[[1]], y = list_samples[[2]])
nonzero <- samples.merged[["RNA"]]@counts > 0
# Keep genes that are expressed in at least 5 cells
tokeep <- rowSums(nonzero) >= 5
genes <- names(tokeep)[tokeep]
samples.merged <- subset(samples.merged, features = genes)

# Perform integration with rpca
list_samples <- SplitObject(samples.merged, split.by = "orig.ident")
list_samples <- lapply(list_samples, SCTransform)

# Select integration features
samples.features <- SelectIntegrationFeatures(object.list = list_samples, nfeatures = 3000)

list_samples <- PrepSCTIntegration(
  object.list = list_samples, anchor.features = samples.features,
  verbose = FALSE
)

# Run PCA before integration
list_samples <- lapply(list_samples, function(x) {
  x <- RunPCA(x, features = samples.features)
})

# Integration considering reciprocal PCA and the SCT normalization
samples.anchors <- FindIntegrationAnchors(
  object.list = list_samples, normalization.method = "SCT",
  anchor.features = samples.features, reduction = "rpca",
  verbose = FALSE
)

samples.integrated <- IntegrateData(
  anchorset = samples.anchors, normalization.method = "SCT",
  verbose = FALSE
)

rm("list_samples")

samples.integrated <- RunPCA(samples.integrated, verbose = TRUE)

# Also check different statistics for the PCs
pct <- samples.integrated[["pca"]]@stdev / sum(samples.integrated[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co1 # 40

# last point where change of % of variation is more than 0.1%.
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
min(co1, co2) # 19

DefaultAssay(samples.integrated) <- "integrated"

# Do it for 20PCs
samples.integrated <- RunUMAP(samples.integrated, dims = 1:20)
samples.integrated <- RunTSNE(samples.integrated, dims = 1:20)
samples.integrated <- FindNeighbors(samples.integrated, dims = 1:20) %>% FindClusters(resolution = 0.7)

Idents(samples.integrated) <- "integrated_snn_res.0.7"

# Separate cluster 7 into Paneth and Goblet
samples.integrated <- FindSubCluster(samples.integrated,
  graph.name = "integrated_snn", cluster = "6", resolution = 0.6
)

samples.integrated$sub.cluster <-
  factor(samples.integrated$sub.cluster,
    levels = c(
      "0", "1", "2", "3", "4", "5", "6_0",
      "6_1", "6_2", "6_3", "7", "8", "9", "10", "11", "12", "13", "14", "15"
    )
  )

Idents(samples.integrated) <- "sub.cluster"

DefaultAssay(samples.integrated) <- "RNA"
samples.integrated <- NormalizeData(samples.integrated)

# Get markers for each of the clusters
annotations <- read_tsv("data/annotation_Mus_v98.tsv")

temp <- FindAllMarkers(samples.integrated, min.pct = 0, logfc.threshold = 0, only.pos = T) %>%
  mutate(specific = pct.1 - pct.2) %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))

temp <- temp %>% filter(p_val_adj < 0.05 & avg_log2FC > 0.25)

write_tsv(temp, file = "Markers_Clusters_subcluster.xls")

ConservedMarkers_list <- vector(mode = "list", length = length(unique(samples.integrated$sub.cluster)))
names(ConservedMarkers_list) <- levels(samples.integrated)
for (i in levels(samples.integrated)) {
  ConservedMarkers_list[[i]] <- FindConservedMarkers(samples.integrated, ident.1 = i, grouping.var = "orig.ident", only.pos = T, min.pct = 0.1, logfc.threshold = 0.1) %>%
    rownames_to_column(var = "gene") %>%
    mutate(
      avg_log2FC = (RANKL_Organoid_avg_log2FC + WT_Organoid_avg_log2FC) / 2,
      specific = ((RANKL_Organoid_pct.1 - RANKL_Organoid_pct.2) + (WT_Organoid_pct.1 - WT_Organoid_pct.2)) / 2,
      cluster = i
    ) %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))
}

ConservedMarkers <- do.call(rbind, ConservedMarkers_list)
write_tsv(ConservedMarkers, file = "ConservedMarkers_subclusters.xls")

# Transfer annotation from reference without distal enterocytes.
reference <- readRDS("data/Haber_2017/Haber_integrated_noDistal.rds")

reference.anchors <- FindTransferAnchors(
  reference = reference, query = samples.integrated, dims = 1:30,
  reference.assay = "integrated",
  query.assay = "integrated",
  normalization.method = "SCT"
)

predictions_simplifiedTA <- TransferData(
  anchorset = reference.anchors,
  refdata = reference$Cluster_SimplifiedTA, dims = 1:30
)

samples.integrated$Clusters_Haber_SimplifiedTA <- predictions_simplifiedTA$predicted.id
samples.integrated$Clusters_Haber_Score_SimplifiedTA <- predictions_simplifiedTA$prediction.score.max

#Only 4 cells are Enterocyte.Progenitor, so we consider them as early looking at annotation
samples.integrated$Clusters_Haber_SimplifiedTA <- factor(samples.integrated$Clusters_Haber_SimplifiedTA,
  levels = c(
    "Enterocyte.Immature.Proximal", "Enterocyte.Mature.Proximal",
    "Enterocyte.Progenitor", "Enterocyte.Progenitor.Early",
    "Enterocyte.Progenitor.Late", "Enteroendocrine", "Goblet",
    "Paneth", "Stem", "TA", "Tuft"
  ),
  labels = c(
    "Enterocyte.Immature.Proximal", "Enterocyte.Mature.Proximal",
    "Enterocyte.Progenitor.Early", "Enterocyte.Progenitor.Early",
    "Enterocyte.Progenitor.Late", "Enteroendocrine", "Goblet",
    "Paneth", "Stem", "TA", "Tuft"
  )
)

#Convert those cells from cluster 7 into microfold, as that celltype was not in the reference
samples.integrated$Clusters_Haber_SimplifiedTA_Microfold <- ifelse(samples.integrated$sub.cluster == "7", "Microfold.cells",
  as.character(samples.integrated$Clusters_Haber_SimplifiedTA)
)

p1 <- DimPlot(samples.integrated)
ggratio <- diff(range(p1$data[, 2])) / diff(range(p1$data[, 1]))

pdf("DimPlot_20PCs.pdf", width = 10)
DimPlot(samples.integrated, group.by = "orig.ident") & theme(legend.position = "top") &
guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 2))) &
coord_fixed(ratio = ggratio)
DimPlot(samples.integrated, group.by = "orig.ident", split.by = "orig.ident", ncol = 2) &
theme(legend.position = "top") & guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 2))) &
coord_fixed(ratio = ggratio)
DimPlot(samples.integrated, group.by = "Phase") & theme(legend.position = "top") &
guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 2))) & 
coord_fixed(ratio = ggratio)
DimPlot(samples.integrated, group.by = "integrated_snn_res.0.7", label = T) & theme(legend.position = "top") &
guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 2))) & coord_fixed(ratio = ggratio)
DimPlot(samples.integrated, group.by = "Clusters_Haber_SimplifiedTA", label = T) & coord_fixed(ratio = ggratio)
DimPlot(samples.integrated, group.by = "Clusters_Haber_SimplifiedTA", split.by = "integrated_snn_res.0.7", ncol = 5) & theme(legend.position = "top") & coord_fixed(ratio = ggratio) & NoLegend()

DefaultAssay(samples.integrated) <- "RNA"
samples.integrated <- NormalizeData(samples.integrated)
FeaturePlot(samples.integrated, features = c("Fabp1", "Agr2", "Tff3", "Defa17", "Slc12a2", "Chgb", "Trpm5", "Lgr5"), order = T, ncol = 4, cols = c("lightgray", "red2")) & coord_fixed(ratio = ggratio)
FeaturePlot(samples.integrated, features = "Clusters_Haber_Score_Simplified", order = T, cols = c("lightgray", "red2")) & coord_fixed(ratio = ggratio)
VlnPlot(samples.integrated, group.by = "integrated_snn_res.0.7", features = c("Fabp1", "Agr2", "Tff3", "Defa17", "Slc12a2", "Chgb", "Trpm5", "Lgr5"), stack = T)
VlnPlot(samples.integrated, group.by = "Clusters_Haber_SimplifiedTA", features = c("Fabp1", "Agr2", "Tff3", "Defa17", "Slc12a2", "Chgb", "Trpm5", "Lgr5"), stack = T)

samples.integrated@meta.data %>% ggplot(aes(y = nFeature_RNA, fill = integrated_snn_res.0.7)) +
  geom_boxplot() +
  theme_classic() +
  scale_y_log10()
samples.integrated@meta.data %>% ggplot(aes(y = percent.mt, fill = integrated_snn_res.0.7)) +
  geom_boxplot() +
  theme_classic()
dev.off()

# Annotate each of the subclusters
samples.integrated$CellType <- factor(samples.integrated$sub.cluster,
  levels = c(
    "0", "1", "2", "3", "4", "5", "6_0", "6_1",
    "6_2", "6_3", "7", "8", "9", "10", "11", "12", "13", "14", "15"
  ),
  labels = c(
    "TA cells", "Stem cells", "TA cells", "Enterocyte Progenitor Late",
    "Enterocyte Mature Proximal", "Enterocyte Mature Proximal", "Goblet cells",
    "Paneth cells", "Goblet cells", "Goblet cells", "Microfold cells", "TA cells",
    "Enterocyte Progenitor Early", "Enterocyte Progenitor Late", "TA cells",
    "Enterocyte Progenitor Late", "Enterocyte Mature Proximal",
    "Enteroendocrine cells", "Tuft cells"
  )
)

samples.integrated$orig.ident <-
  factor(samples.integrated$orig.ident,
    levels = c("WT_Organoid", "RANKL_Organoid"),
    labels = c("Ctrl", "RANKL")
  )

# Control organoids have no microfold cells
# thus need to convert those 6 cells in cluster 6 to enterocytes
samples.integrated$CellType <- ifelse(samples.integrated$sub.cluster == "7" &
  samples.integrated$orig.ident == "Ctrl",
"Enterocyte Mature Proximal", as.character(samples.integrated$CellType)
)

samples.integrated$FinalCTCondition <- paste(samples.integrated$CellType,
  samples.integrated$orig.ident,
  sep = "_"
)

# Reorder celltype levels as desired
celltypes <- c(
  "Stem cells", "TA cells", "Enterocyte Progenitor Early",
  "Enterocyte Progenitor Late",
  "Enterocyte Mature Proximal", "Microfold cells",
  "Enteroendocrine cells", "Tuft cells", "Paneth cells", "Goblet cells"
)

samples.integrated$CellType <- factor(samples.integrated$CellType,
  levels = celltypes
)

write.table(samples.integrated@meta.data, file = "metadata.xls", row.names = T)


saveRDS(samples.integrated, file = "Organoids.rds")
