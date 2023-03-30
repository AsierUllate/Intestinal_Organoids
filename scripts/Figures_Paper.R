library(Seurat)
library(tidyverse)
library(future)
library(pheatmap)
library(RColorBrewer)
library(msigdbr)
library(scales)
library(UCell)
library(ggprism)
options(future.globals.maxSize = 12000 * 1024^2)
options(bitmapType = "cairo")

plan("multiprocess", workers = 6)


# Load integrated sample
samples.integrated <- readRDS("Organoids.rds")
DefaultAssay(samples.integrated) <- "RNA"
samples.integrated <- NormalizeData(samples.integrated)

colors_samples <- c("Ctrl" = "#0000CC", "RANKL" = "#D00004")

#Get ratio for representing UMAPs
p1 <- DimPlot(samples.integrated)
ggratio <- diff(range(p1$data[, 2])) / diff(range(p1$data[, 1]))

#Get DimPlot colors for simplified TA in Haber
col_celltypes <- c(RColorBrewer::brewer.pal(n = 10, name = "Paired"))

names(col_celltypes) <- c("Stem cells", "TA cells", "Enterocyte Progenitor Early",
"Enterocyte Progenitor Late",
"Microfold cells", "Enterocyte Mature Proximal",
"Enteroendocrine cells", "Tuft cells", "Paneth cells", "Goblet cells")


# Generate DimPlot with the new annotation

pdf("DimPlot_CellType_ratio.pdf", width = 12)
DimPlot(samples.integrated, group.by = "CellType", label = F,
label.size = 5, cols = col_celltypes, repel = T) +
ggtitle("Epithelial cell types based on unsupervised clustering") +
theme(plot.title = element_text(hjust = 0.5))  + coord_fixed(ratio = ggratio)
dev.off()

##Extended data Figure 5 Panel A
pdf("DimPlot_CellType_SplitSample.pdf", width = 18)
DimPlot(samples.integrated, group.by = "CellType", split.by = "orig.ident", label = F,
cols = col_celltypes, repel = F)  & coord_fixed(ratio = ggratio) &
theme(plot.title = element_blank()) 
dev.off()

pdf("DimPlot_CellType_SplitSample.pdf", width = 10)
DimPlot(samples.integrated, group.by = "CellType", split.by = "orig.ident", label = F,
cols = col_celltypes, repel = F)  & coord_fixed(ratio = ggratio) &
theme(plot.title = element_blank(), legend.position = "bottom",
legend.justification = "center", legend.text = element_text(size = 16))  &
guides(color = guide_legend(nrow = 5, byrow = FALSE, override.aes = list(size = 5)))
dev.off()

#ViolinPlots Figure 1e
genes <- c("Birc2", "Birc3", "Tnfaip3", "Bcl2", "Bcl2l1")

pdf("VlnPlot_Bircs.pdf", width = 7)
VlnPlot(samples.integrated, features = genes,
group.by = "CellType", split.by = "orig.ident",
stack = T, flip = T) + scale_fill_manual(values = colors_samples) +
theme(axis.text.x = element_text(angle =  45)) + xlab(element_blank()) + 
ylab("Log-normalized Expression")
dev.off()

#ViolinPlots Figure 1g
bmp_genes <- c("Bmp2", "Id2", "Id3")
pdf("VlnPlot_Bmps.pdf", width = 7)
VlnPlot(samples.integrated, features = bmp_genes,
group.by = "CellType", split.by = "orig.ident", stack = T, flip = T) +
scale_fill_manual(values = colors_samples) +
theme(axis.text.x = element_text(angle =  45)) + xlab(element_blank()) + 
ylab("Log-normalized Expression")
dev.off()


#VlnPlot from Rank in Control Sample for Extended Figure 5c
pdf("VlnPlot_Tnfrsf11a_Ctrl.pdf", width = 7)
VlnPlot(subset(samples.integrated, subset = orig.ident == "Ctrl"),
features = "Tnfrsf11a", cols = col_celltypes, group.by = "CellType",
pt.size = 0.1) + theme(axis.text.x = element_text(angle =  45)) +
xlab(element_blank()) + 
ylab("Rank log-normalized Expression") + NoLegend()
dev.off()

#Heatmap from markers Yum and Biton markers of differenti celltypes
# Extended Data Figure 5b
markers <- read_tsv(
"data/Markers_Yum_Biton_epithelial.tsv", col_names = F)

markers_list <- split(markers$X1, markers$X2)
DefaultAssay(samples.integrated) <- "RNA"
samples.integrated <- ScaleData(samples.integrated, features = markers$X1)

#Get data for pheatmap
i <- "CellType"
genes <- FetchData(samples.integrated,
vars = c(markers$X1, all_of(i), "orig.ident"),
slot = "scale.data") %>% rownames_to_column(var = "cellname")

genes <- genes %>% arrange_at(i)
annotations <- data.frame(row.names = genes$cellname,
Celltype = genes[, i],
Sample = genes[, "orig.ident"])

#annotation of genes
genes_short <- genes %>% dplyr::select(!c(all_of(i), "orig.ident")) %>%
column_to_rownames(var = "cellname")
annotation_genes <- data.frame(row.names = colnames(genes_short),
Markers = markers$X2[match(colnames(genes_short), markers$X1)])

#annotation colors
column_cols <- col_celltypes

Status_cols <- hue_pal()(length(unique(annotation_genes$Markers)))
names(Status_cols) <- unique(annotation_genes$Markers)
Sample_cols <- colors_samples
ann_colors  <- list(Celltype  = column_cols,
                    Markers = Status_cols,
                    Sample = Sample_cols)
genes_short <- t(MinMax(genes_short, min = -3, max = 3))

pdf("Heatmaps_markers_Celltype.pdf", width = 10,
height = 12)

pheatmap(genes_short, show_colnames = F, cluster_cols = F,
color = colorRampPalette(rev(brewer.pal(n = 9,name = "RdBu")))(50),
cluster_rows = F, show_rownames = T,
annotation_col = annotations, annotation_row = annotation_genes,
annotation_colors = ann_colors,
annotation_names_row = F, main = i,
gaps_row = as.numeric(cumsum(table(annotation_genes$Markers)[
  unique(annotation_genes$Markers)])),
  gaps_col = as.numeric(cumsum(table(annotations$Celltype))))
dev.off()

######
# Check different signatures for Extended Data Fig.6
#####

#Load the gene sets to compare
load(file = "data/gene_sets.Rdata")

#Try with UCell
samples.integrated <- AddModuleScore_UCell(samples.integrated,
features = gene_sets, assay = "RNA")

signature_names_UCell <- paste0(names(gene_sets), "_UCell")

info_UCell <- samples.integrated@meta.data[, signature_names_UCell]

save(info_UCell, signature_names_UCell,
file = "data/UCell.RData")

load(file.path(project_dir, "data/UCell.RData"))


#####
# Do it for all wanted columns from the meta.data, it can be split in several folders if wanted
####
comp <- "UCell"
if (!dir.exists(file.path(comp))) {
dir.create(file.path(comp))
}
# Generate all the folders the first time
if (!dir.exists(file.path(comp, "VlnPlot_GeneSets_significance_adj"))) {
dir.create(file.path(comp, "VlnPlot_GeneSets_significance_adj"))
}
#Get all signature names in a dataframe
df <- FetchData(samples.integrated, vars = c(signature_names_UCell,
"orig.ident", "CellType"))
#df <- df %>% filter(CellType != "Microfold cells")
df$CellType <- droplevels(df$CellType)
df$orig.ident <- factor(df$orig.ident)
#In long format to filter and get those that we need to
df_long <- df %>% pivot_longer(cols = all_of(signature_names_UCell), names_to = "geneset")

names(signature_names_UCell) <- c("BMP pathway", "Anti-apoptotic pathway", "NF-kB pathway", "ERK/MAPK targets")
#Generate plots for each of the gene sets
# Put only significance
for (sig in names(signature_names_UCell)) {
    ii <- signature_names_UCell[sig]
    # Plot adjusted p.values
    # Correct with BH

    #Get dataframe with pvalues for the gene set of interest
    #Add some data for Ctrl Microfold cells
    df_long_tmp <- as.data.frame(t(c("Ctrl", "Microfold cells", ii, 0.0)))
    colnames(df_long_tmp) <- colnames(df_long)
    df_long_tmp$value <- as.numeric(df_long_tmp$value)

    df_p_val <- df_long %>%
    filter(geneset == ii) %>% rbind(., df_long_tmp) %>%
    rstatix::group_by(CellType) %>%
    rstatix::wilcox_test(value ~ orig.ident) %>%
    rstatix::adjust_pvalue(p.col = "p", method = "BH") %>%
    rstatix::add_significance(p.col = "p.adj") %>%
    rstatix::add_xy_position(x = "CellType", dodge = 0.8)
    #Reduce to 2 digits
    df_p_val$p.adj <- signif(df_p_val$p.adj, digits = 2)
    df_p_val$p.adj.signif[df_p_val$CellType == "Microfold cells"] <- NA

    #Get Violin Plot with geom_violin as if not it doesn't work with add_pvalue
    p1 <- samples.integrated@meta.data %>%
    ggplot(aes_string(x = "CellType",
    y = ii)) + geom_violin(aes(fill = orig.ident), scale = "width") +
    scale_fill_manual(values = colors_samples) +
    #geom_jitter position is as in Seurat
    geom_jitter(aes(fill = orig.ident), size = 0.05,
    position = position_jitterdodge(jitter.width = 0.4,
    dodge.width = 0.9), show.legend = FALSE) + theme_classic() +
    theme(axis.text.x = element_text(angle =  45, size = 18, hjust = 1),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title.y = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_blank()) + ylab("UCell score") +
    xlab(element_blank()) + ggtitle(sig)

    # Position in highest y.position in violins
    df_p_val$y.position <- max(df_p_val$y.position)
    pdf(file.path(comp,
    sprintf("VlnPlot_GeneSets_significance_adj/VlnPlot_%s_pval_adj_high.pdf", ii)), height = 8, width = 9)
    p2 <- p1 + add_pvalue(df_p_val, xmin = "xmin", label.size = 7,
    xmax = "xmax", label = "p.adj.signif", tip.length = 0, bracket.size = 0, bracket.colour = "white")
    print(p2)
    dev.off()
}


#####
#Get stastistical p-values from Wilcoxon test
####
list_pvals <- vector(mode = "list", length = length(signature_names_UCell))
names(list_pvals) <- signature_names_UCell
list_pval_adj <- list_pvals
for (i in names(list_pvals)) {
cluster <- c(unique(as.vector(as.matrix(df[, "CellType"]))))
cluster <- sort(cluster)
cluster <- cluster[cluster != "Microfold cells"]
pval <- vector(mode = "double", length = length(cluster))
names(pval) <- cluster
for (lev in levels(factor(df[, "CellType"]))) {
    if (lev != "Microfold cells"){
        pval[lev] <- wilcox.test(get(i) ~ orig.ident,
        data = df, subset = CellType == lev)$p.value
    }
}
list_pvals[[i]] <- data.frame(cluster = cluster,
pvalue = pval, signature = i)
list_pval_adj[[i]] <- data.frame(cluster = cluster,
pvalue = p.adjust(pval, method = "BH"), signature = i)
}
df_pvals <- do.call(rbind, list_pvals)
df_pvals_adj <- do.call(rbind, list_pval_adj)
write_tsv(df_pvals, file =
file.path(comp, "pvalues_wilcoxtest.tsv"))
write_tsv(df_pvals_adj, file =
file.path(comp, "pvalues_wilcoxtest_BHadj.tsv"))

df_pvals_wide <- pivot_wider(df_pvals,
names_from = signature, values_from = pvalue)
write_tsv(df_pvals_wide,
file = file.path(comp, "pvalues_wilcoxtest_wide.tsv"))

df_pvals_adj_wide <- pivot_wider(df_pvals_adj,
names_from = signature, values_from = pvalue)
write_tsv(df_pvals_adj_wide,
file = file.path(comp, "pvalues_wilcoxtest_BHadj_wide.tsv"))


Intestine_Paper <- readRDS("data/Fujii_2018/Fujii.rds")

cols_clusters <- RColorBrewer::brewer.pal("Paired", n = 12)

#Extended Data Fig. 16a
pdf("UMAP_CellType_Intestine_25PCs_epithelialonly_UMAP.pdf")
DimPlot(Intestine_Paper, label = F, repel = F, cols = c(cols_clusters),
group.by = "CellType") & coord_fixed(ratio = ggratio)
dev.off()

#Extended Data Fig. 16b
pdf("VlnPlot_RANK_Celltype_epithelialonly.pdf", width = 10)
p1 <- VlnPlot(Intestine_Paper, features = "TNFRSF11A", cols = c(cols_clusters),
group.by = "CellType", fill.by = "ident", pt.size = 0.1) &
theme(axis.text.x = element_text(angle =  45)) & xlab(element_blank()) &
ylab("Log-normalized Expression") & ggtitle("RANK") & NoLegend()

p2 <- VlnPlot(Intestine_Paper, features = "TNFSF11", cols = c(cols_clusters),
group.by = "CellType", fill.by = "ident", pt.size = 0.1) &
theme(axis.text.x = element_text(angle =  45)) & xlab(element_blank()) &
ylab("Log-normalized Expression") & ggtitle("RANKL") & NoLegend()
p1 | p2
dev.off()