library(Seurat)
library(dplyr)
library(reticulate)
library(cowplot)
library(ggplot2)
library(viridis)
library(tidyr)
library(magrittr)
library(reshape2)
library(readr)
library(stringr)
library(RColorBrewer)
library(pheatmap)
library(cols4all)
library(ComplexHeatmap)
library(circlize)
library(leidenbase)
library(patchwork)
library(monocle3)
library(pROC)
library(ggrepel)

results_dir <- "results"
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

### reading 10X data obtained from GEO database (GSE198204) ====================
p1.AF <- Read10X("GSE198204_human/AF1/")
p1.AF <- CreateSeuratObject(counts = p1.AF, project = "Patient1.AF", min.features = 200)
p1.AF$percent.mt <- PercentageFeatureSet(p1.AF, pattern = "^MT-")
p1.AF$disease <- "AF"
p1.nor <- Read10X("GSE198204_human/No AF1/")
p1.nor <- CreateSeuratObject(counts = p1.nor, project = "Patient1.normal", min.features = 200)
p1.nor$percent.mt <- PercentageFeatureSet(p1.nor, pattern = "^MT-")
p1.nor$disease <- "Normal"
p2.AF <- Read10X("GSE198204_human/AF2/")
p2.AF <- CreateSeuratObject(counts = p2.AF, project = "Patient2.AF", min.features = 200)
p2.AF$percent.mt <- PercentageFeatureSet(p2.AF, pattern = "^MT-")
p2.AF$disease <- "AF"
p2.nor <- Read10X("GSE198204_human/No AF2/")
p2.nor <- CreateSeuratObject(counts = p2.nor, project = "Patient2.normal", min.features = 200)
p2.nor$percent.mt <- PercentageFeatureSet(p2.nor, pattern = "^MT-")
p2.nor$disease <- "Normal"
p3.AF <- Read10X("GSE198204_human/AF3/")
p3.AF <- CreateSeuratObject(counts = p3.AF, project = "Patient3.AF", min.features = 200)
p3.AF$percent.mt <- PercentageFeatureSet(p3.AF, pattern = "^MT-")
p3.AF$disease <- "AF"
p3.nor <- Read10X("GSE198204_human/No AF3/")
p3.nor <- CreateSeuratObject(counts = p3.nor, project = "Patient3.normal", min.features = 200)
p3.nor$percent.mt <- PercentageFeatureSet(p3.nor, pattern = "^MT-")
p3.nor$disease <- "Normal"

scRNAlist <- list(p1.AF, p1.nor, p2.AF, p2.nor, p3.AF, p3.nor)
scRNA_harmony <- merge(scRNAlist[[1]], y = scRNAlist[-1])

## Figure 4A
VlnPlot(scRNA_harmony, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

### data preprocessing using SCTransform
pbmc <- scRNA_harmony
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 &
                          nCount_RNA > 500 & percent.mt < 15)
pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
ElbowPlot(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.1)
pbmc <- RunUMAP(pbmc, dims = 1:20)

## Figure 4B
DimPlot(pbmc, label = TRUE)
saveRDS(pbmc, file = "pbmc.rds")

### cell type annotation
seurat <- readRDS("pbmc.rds")
seurat <- PrepSCTFindMarkers(seurat)
markers <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
seurat$cell_type <- recode(as.character(seurat$seurat_clusters),
    "0"  = "CD4+ T cells",
    "1"  = "NK cells",
    "2"  = "Macrophages",
    "3"  = "CD8+ T cells",
    "4"  = "Inflammatory monocytes",
    "5"  = "Endothelial cells",
    "6"  = "Fibroblasts",
    "7"  = "Vascular smooth muscle cells",
    "8"  = "Endothelial cells(activated)",
    "9"  = "Neutrophils",
    "10" = "B cells"
)

## Figure 4C
DimPlot(seurat, group.by = "cell_type", label = TRUE, pt.size = 0.5)
seurat <- JoinLayers(seurat)
seurat <- NormalizeData(seurat, assay = "RNA")
seurat <- FindVariableFeatures(seurat, assay = "RNA")
seurat <- ScaleData(seurat, assay = "RNA")
saveRDS(seurat, file = "seurat.rds")

## Figure 4D
marker <- FindAllMarkers(
  object = seurat,
  test.use = "wilcox",
  only.pos = TRUE,
  logfc.threshold = 0.5,
  min.pct = 0.5
)
all.markers <- subset(marker, p_val < 0.05)
top5 <- all.markers |>
  group_by(cluster) |>
  slice_max(n = 3, order_by = avg_log2FC) |>
  ungroup()
genes <- unique(top5$gene)
aver_dt <- AverageExpression(seurat, features = genes,
                             group.by = "cell_type", slot = "data")$RNA
aver_dt_ordered <- aver_dt[top5$gene, ]
aver_dtt <- t(scale(t(aver_dt_ordered)))

gene_anno <- data.frame(cluster = top5$cluster, row.names = top5$gene)
celltype_colors <- c(
  "Macrophages"                  = "#FFB3BA",
  "T cells"                      = "#BAE1FF",
  "Endothelial cells"            = "#BAFFC9",
  "Fibroblasts"                  = "#FFFFBA",
  "Atrial cardiomyocytes"        = "#FFB7CE",
  "Vascular smooth muscle cells" = "#E0BBFF",
  "Neutrophils"                  = "#BFFCC6",
  "Swwann cells"                 = "#FFDFBA",
  "Lymphatic endothelial cells"  = "#A8D8EA"
)
cluster_levels <- unique(top5$cluster)
cluster_palette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00",
                     "#984EA3", "#A65628", "#F781BF", "#999999", "#66C2A5")
cluster_colors <- setNames(cluster_palette[seq_along(cluster_levels)], cluster_levels)
mycol <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#D73027"))
col_anno <- HeatmapAnnotation(
  cell_type = colnames(aver_dtt),
  col = list(cell_type = celltype_colors),
  show_annotation_name = FALSE,
  gp = gpar(col = "white", lwd = 2)
)
row_anno <- rowAnnotation(
  cluster = gene_anno[rownames(aver_dtt), "cluster"],
  col = list(cluster = cluster_colors),
  show_annotation_name = FALSE,
  gp = gpar(col = "white", lwd = 1.5)
)
Heatmap(
  aver_dtt,
  name = "Scaled\nexpression",
  col = mycol,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  column_names_side = "top",
  column_names_rot = 60,
  row_names_gp = gpar(fontsize = 10, fontface = "italic"),
  rect_gp = gpar(col = "white", lwd = 1.5),
  top_annotation = col_anno,
  left_annotation = row_anno,
  row_title = NULL,
  border = TRUE
)

### obtaining DEGs in Macrophages between AF and Normal ========================
Idents(seurat) <- seurat$disease
DefaultAssay(seurat) <- "RNA"
Macro <- subset(seurat, cell_type %in% c("Macrophages", "Inflammatory monocytes"))

markers_af_vs_normal <- FindMarkers(Macro,
                                    ident.1 = "AF",
                                    ident.2 = "Normal",
                                    min.pct = 0.25,
                                    logfc.threshold = 0.25)
markers_af_vs_normal$gene <- rownames(markers_af_vs_normal)
log2FC_cutoff <- 0.25
pval_cutoff  <- 0.05
markers_af_vs_normal$sig <- "NS"
markers_af_vs_normal$sig[markers_af_vs_normal$avg_log2FC >  log2FC_cutoff &
                           markers_af_vs_normal$p_val_adj < pval_cutoff] <- "Up"
markers_af_vs_normal$sig[markers_af_vs_normal$avg_log2FC < -log2FC_cutoff &
                           markers_af_vs_normal$p_val_adj < pval_cutoff] <- "Down"
rps26_df <- markers_af_vs_normal[markers_af_vs_normal$gene == "RPS26", , drop = FALSE]

pbmc <- Macro
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 &
                         nCount_RNA > 500 & percent.mt < 15)
pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.1)
pbmc <- RunUMAP(pbmc, dims = 1:20)
Macro <- pbmc

## Figure 4E
volcano_plot <- ggplot(markers_af_vs_normal,
                       aes(x = avg_log2FC, y = -log10(p_val_adj), color = sig)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c(Up = "#D73027", Down = "#2166AC", NS = "grey80")) +
  geom_vline(xintercept = c(-log2FC_cutoff, log2FC_cutoff),
             linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(pval_cutoff),
             linetype = "dashed", color = "grey50") +
  labs(x = "log2 fold change", y = "-log10(adjusted p-value)", color = NULL) +
  theme_classic(base_size = 14) +
  theme(panel.grid = element_blank(), legend.position = "top")

volcano_plot_RPS26 <- volcano_plot +
  geom_point(data = rps26_df, shape = 1, color = "black", size = 3.5, stroke = 1) +
  geom_label_repel(data = rps26_df, aes(label = gene),
                   size = 3, max.overlaps = 20, segment.color = "grey50",
                   min.segment.length = 0, box.padding = 0.5, seed = 42,
                   label.size = 0.4, label.padding = unit(0.25, "lines"),
                   label.r = unit(0.15, "lines"), fill = "white") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14))
volcano_plot_RPS26

