library(Seurat)
library(dplyr)
library(sceasy)
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
###devtools::install_github("cellgeni/schard")
library(schard)
###1.h5ad data reading
data <- h5ad2seurat("GSE224959_atrial fibrillation_cellxgene.h5ad")
data$percent.mt <- PercentageFeatureSet(data , pattern = "^MT-")
scRNAlist <- list(data=data)
pbmc <- merge(scRNAlist[[1]], y=c(scRNAlist[[2]]))
for(i in seq(4,length(scRNAlist))){
  pbmc <- merge(pbmc, y=c(scRNAlist[[i]]))
}

### 2. data preprocessing through SCTransform
pbmc <- subset(pbmc,subset=nFeature_RNA > 200 & nFeature_RNA < 2500 & nCount_RNA > 500 & percent.mt < 15)
 ## Figure 1A
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#options(future.globals.maxSize = 16 * 1024^3) 
pbmc <- SCTransform(pbmc,vars.to.regress="percent.mt",verbose=FALSE)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
ElbowPlot(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.1)
pbmc <- RunUMAP(pbmc, dims = 1:20)
DimPlot(pbmc, label = TRUE)
saveRDS(pbmc, file = "pbmc.rds")

### 3. cell type annotation
seurat <- readRDS("pbmc.rds")
seurat <- PrepSCTFindMarkers(seurat)
markers <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
seurat$cell_type <- recode(as.character(seurat$seurat_clusters),
    "0" = "CD4+ T cells",
    "1" = "NK cells",
    "2" = "Macrophages",
    "3" = "CD8+ T cells",
    "4" = "Inflammatory monocytes",
    "5" = "Endothelial cells",
    "6" = "Fibroblasts",
    "7" = "Vascular smooth muscle cells",
    "8" = "Endothelial cells(activated)",
    "9" = "Neutrophils",
    "10" = "B cells"
)
## Figure 1B
DimPlot(seurat, group.by = "cell_type", label = TRUE, pt.size = 0.5)

seurat <- JoinLayers(seurat)
seurat <- NormalizeData(seurat, assay = "RNA")
seurat <- FindVariableFeatures(seurat, assay = "RNA")
seurat <- ScaleData(seurat, assay = "RNA")
saveRDS(seurat, file = "seurat.rds")

## Figure 1C
seurat <- readRDS("seurat.rds")
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
aver_dt <- AverageExpression(
  seurat,
  features = genes,
  group.by = "cell_type",  
  slot = "data"
)$RNA

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
### 4. cell type composition analysis
compare_df <- prop_df |>
  select(disease, cell_type, n, prop) |>
  pivot_wider(
    id_cols = cell_type,
    names_from = disease,
    values_from = c(n, prop)
  ) |>
  mutate(
    diff = `prop_AF/MVI` - `prop_normal`,
    across(where(is.numeric), ~ round(.x, 2))
  ) |>
  arrange(desc(diff))

compare_df |>
  mutate(
    `prop_normal` = sprintf("%.2f%%", `prop_normal`),
    `prop_AF/MVI` = sprintf("%.2f%%", `prop_AF/MVI`)
  ) |>
  print(n = Inf)
prop_df$disease <- factor(prop_df$disease,
  levels = c("Normal", "atrial fibrillation || mitral valve insufficiency"),
  labels = c("Normal", "AF/MVI")
)
ggplot(prop_df, aes(x = disease, y = prop, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.6, color = "white", linewidth = 0.3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02)), limits = c(0, 100)) +
  scale_fill_brewer(palette = "Set3") +
  labs(
    x = NULL,
    y = "Cell proportion (%)",
    fill = "Cell type",
    title = "Cell type composition by disease"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "right"
  )


#### 5. Differential expression analysis in Macrophages subcluster
Idents(seurat) <- seurat$disease
DefaultAssay(seurat) <- "RNA"
Macro <- subset(seurat, cell_type == c("Macrophages")|cell_type == "Inflammatory monocytes")
markers_af_vs_normal <- FindMarkers(Macro, 
                                    ident.1 = "atrial fibrillation || mitral valve insufficiency", 
                                    ident.2 = "normal",
                                    min.pct = 0.25,
                                    logfc.threshold = 0.25)
markers_af_vs_normal$gene <- rownames(markers_af_vs_normal)
pbmc <- Macro
pbmc <- subset(pbmc,subset=nFeature_RNA > 200 & nFeature_RNA < 2500 & nCount_RNA > 500 & percent.mt < 15)
#options(future.globals.maxSize = 16 * 1024^3)
#pbmc <- SCTransform(pbmc,vars.to.regress="percent.mt",verbose=FALSE)
pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt",method = "glmGamPoi",
                     return.only.var.genes = TRUE, verbose = FALSE)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
ElbowPlot(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.1)
pbmc <- RunUMAP(pbmc, dims = 1:20)
## Figure 1E
DimPlot(pbmc, label = TRUE)
Macro <- pbmc
saveRDS(Macro, file = "Macro.rds")
## Figure 2A
log2FC_cutoff <- 0.25
pval_cutoff  <- 0.05
markers_af_vs_normal <- markers_af_vs_normal |>
  mutate(
    gene = rownames(markers_af_vs_normal),
    sig = case_when(
      avg_log2FC >  log2FC_cutoff & p_val_adj < pval_cutoff ~ "Up",
      avg_log2FC < -log2FC_cutoff & p_val_adj < pval_cutoff ~ "Down",
      TRUE ~ "NS"
    )
  )
top15 <- markers_af_vs_normal |>
  mutate(score = -log10(p_val_adj) * abs(avg_log2FC)) |>
  arrange(desc(score)) |>
  slice_head(n = 15)
volcano_plot <- ggplot(markers_af_vs_normal,
                       aes(x = avg_log2FC, y = -log10(p_val_adj), color = sig)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c(Up = "#D73027", Down = "#2166AC", NS = "grey80")) +
  geom_vline(xintercept = c(-log2FC_cutoff, log2FC_cutoff),
             linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(pval_cutoff),
             linetype = "dashed", color = "grey50") +
  labs(
    x = "log2 fold change",
    y = "-log10(adjusted p-value)",
    color = NULL,
    title = "AF/MVI vs Normal (Macrophages)"
  ) +
  theme_classic(base_size = 14) +
  theme(panel.grid = element_blank(), legend.position = "top")
volcano_plot_border <- volcano_plot +
  geom_label_repel(
    data = top15,
    aes(label = gene),
    size = 3,
    max.overlaps = 20,
    segment.color = "grey50",
    min.segment.length = 0,
    box.padding = 0.5,
    seed = 42,
    label.size = 0.4,
    label.padding = unit(0.25, "lines"),
    label.r = unit(0.15, "lines"),
    fill = "white"
  ) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )
volcano_plot_border


### 6. Subcluster function enrichment 
Macro <- readRDS("Macro.rds")
Macro$subcluster <- Idents(Macro)
Macro <- subset(Macro,subcluster %in% c(0,1,2))
library(clusterProfiler)
library(org.Hs.eg.db)
Idents(Macro) <- "subcluster"
macro_markers <- FindAllMarkers(
  Macro,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
sig_markers <- macro_markers |>
  dplyr::filter(p_val_adj < 0.05)
gene_list <- split(sig_markers$gene, sig_markers$cluster)
entrez_list <- lapply(gene_list, function(genes) {
  bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
})
kegg_results <- lapply(entrez_list, function(entrez_ids) {
  enrichKEGG(
    gene = entrez_ids,
    organism = "hsa",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
})
readRDS(Macro,file = "Macro.rds")

### 7. Macrophage pseudotime analysis
Macro <- readRDS("Macro.rds")
data  <-  GetAssayData(Macro,  layer = 'counts')
cell_metadata  <-  Macro@meta.data
gene_annotation  <-  data.frame(gene_short_name = rownames(data))
rownames(gene_annotation)  <-  rownames(data)
cds  <-  new_cell_data_set(data,
                            cell_metadata = cell_metadata,
                            gene_metadata = gene_annotation)
cds  <-  preprocess_cds(cds, num_dim = 100)
plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds,
                        reduction_method = "UMAP",
                        preprocess_method = "Aligned",
                        umap.n_neighbors = 40,
                        umap.min_dist = 0.3)
cds <- cluster_cells(cds, resolution = 0.01)
cds <- learn_graph(cds, learn_graph_control = list(
  euclidean_distance_ratio = 1.8,
  geodesic_distance_ratio = 0.4,
  minimal_branch_len = 30,
  rann.k = 60,
  L1.gamma = 0.05
))
cds = order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=F,
           graph_label_size=1.5,
           group_label_size=4,cell_size=1)
saveRDS(cds,file = "Macro_pseudotime.rds")