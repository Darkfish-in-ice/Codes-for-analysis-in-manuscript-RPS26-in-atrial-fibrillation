library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
library(celldex)
library(SingleR)
library(clustree)
library(harmony)
library(assertthat)
library(SCpubr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(enrichplot)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(scTenifoldKnk)
### 1. RPS26 Virtual Knockout in Macrophages
target_gene <- "RPS26"
Macro <- readRDS("Macro.rds")
count_matrix <- GetAssayData(Macro, layer = "counts")
pbmc <- FindVariableFeatures(Macro, selection.method = "vst", nfeatures = 10000)
high_var_genes <- VariableFeatures(pbmc)[1:1000]
input_data <- as.data.frame(count_matrix[c(high_var_genes,target_gene),])
ko_result <- scTenifoldKnk(
   countMatrix = input_data,
   gKO = target_gene,
   qc = FALSE,                   
   nc_nNet = 10,
   nc_nCells = 500
)
tensor_networks <- ko_result$tensorNetworks
manifold_alignment <- ko_result$manifoldAlignment
diff_regulation <- ko_result$diffRegulation
diff <- diff_regulation$gene
diff_entrez<-bitr(
  diff,
  fromType='SYMBOL',
  toType='ENTREZID',
  OrgDb='org.Hs.eg.db'
)
KEGG_enrich<-clusterProfiler::enrichKEGG(gene=diff_entrez$ENTREZID,
                                         organism="hsa",
                                         pvalueCutoff=0.05,
                                         qvalueCutoff=0.05,
                                         pAdjustMethod="BH",
                                         minGSSize=10,
                                         maxGSSize=500)
KEGG_enrich<-setReadable(KEGG_enrich,
                         OrgDb="org.Hs.eg.db",
                         keyType='ENTREZID')
KEGG_result<-KEGG_enrich@result
write.csv(diff_regulation, file = "diff_regulation.csv")


## Figure 5C Volcano plot
log2fc_cutoff <- 1
pval_cutoff <- 0.05
vol_df <- diff_regulation |>
  mutate(
    log2FC = log2(FC),
    neg_log10p = -log10(pmax(p.adj, 1e-300)),
    sig = case_when(
      p.adj < pval_cutoff & log2FC >  log2fc_cutoff ~ "Down in KO",
      p.adj < pval_cutoff & log2FC < -log2fc_cutoff ~ "Up in KO",
      TRUE ~ "NS"
    )
  )
top15 <- vol_df |>
  mutate(score = neg_log10p * abs(log2FC)) |>
  arrange(desc(score)) |>
  slice_head(n = 15)
volcano_ko <- ggplot(vol_df, aes(x = log2FC, y = neg_log10p, color = sig)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c("Up in KO" = "#2166AC", "Down in KO" = "#D73027", "NS" = "grey80")) +
  scale_y_continuous(0, 150) +
  geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(pval_cutoff), linetype = "dashed", color = "grey50") +
  geom_text_repel(
    data = top_up,
    aes(label = gene),
    size = 3,
    max.overlaps = 20,
    segment.color = "grey50",
    min.segment.length = 0,
    box.padding = 0.5,
    seed = 42
  ) +
  labs(
    x = "log2(FC) (WT/KO)",
    y = "-log10(adjusted p-value)",
    color = NULL,
    title = "RPS26 Virtual Knockout"
  ) +
  theme_classic(base_size = 14) +
  theme(panel.grid = element_blank(), legend.position = "top")

volcano_ko


### 2. RPS26 Virtual Overexpression in Macrophages
library(scTenifoldNet)
library(Seurat)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(MASS)
gOE <- "RPS26"
rna_bind <- readRDS("Macro.rds")
rna_bind[["percent.mt"]] <- PercentageFeatureSet(rna_bind, pattern = "^MT-")
VlnPlot(rna_bind,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  pt.size = 0.1,
  combine = FALSE
)
rna_bind_qc <- subset(
  rna_bind,
  subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5
)
count_matrix <- LayerData(rna_bind_qc, assay = "RNA", layer = "counts")
count_matrix <- count_matrix[rowMeans(count_matrix != 0) >= 0.05, ]
if (!gOE %in% rownames(count_matrix)) {
  stop(sprintf("Target gene '%s' not found.", gOE))
}
set.seed(1234)
n_cores <- parallel::detectCores()
WT_networks <- makeNetworks(
  X = count_matrix,
  nNet = 10,
  nCells = 500,
  nComp = 3,
  nCores = n_cores
)
WT_tensor <- tensorDecomposition(WT_networks, K = 3)
WT_network <- WT_tensor$X
OE_network <- WT_network
OE_network[gOE, ] <- OE_network[gOE, ] * 10
MA_results <- manifoldAlignment(
  t(as.matrix(WT_network)),
  t(as.matrix(OE_network))
)
geneList <- rownames(MA_results)
geneList <- gsub("^X_", "", geneList[grepl("^X_", geneList)])
nGenes <- length(geneList)
distances <- sapply(seq_len(nGenes), function(i) {
  stats::dist(rbind(MA_results[i, ], MA_results[i + nGenes, ]))
})
lambda_val <- try({
  boxcox_fit <- MASS::boxcox(distances[distances > 0] ~ 1, plotit = FALSE)
  boxcox_fit$x[which.max(boxcox_fit$y)]
}, silent = TRUE)
if (inherits(lambda_val, "try-error")) {
  transformed_distances <- distances
} else if (lambda_val < 0) {
  transformed_distances <- 1 / (distances^lambda_val)
} else {
  transformed_distances <- distances^lambda_val
}
Z_scores <- as.numeric(scale(transformed_distances))
mean_bg_dist_sq <- mean(distances[-which(geneList == gOE)]^2)
FC <- as.numeric(distances^2 / mean_bg_dist_sq)
pValues <- pchisq(q = FC, df = 1, lower.tail = FALSE)
pAdjusted <- p.adjust(pValues, method = "fdr")
dr_results <- data.frame(
  gene = geneList,
  distance = as.numeric(distances),
  Z = Z_scores,
  FC = FC,
  p.value = pValues,
  p.adj = pAdjusted
) |>
  arrange(desc(Z))
dr_results_plot <- dr_results |>
  mutate(
    significance = ifelse(p.adj < 0.05, "Significant", "Not Significant"),
    label = if_else(p.adj < 0.05 & Z > 3, gene, NA_character_)
  )

## Figure 5F Volcano plot
diff_reg_plot <- ggplot(dr_results_plot, aes(x = Z, y = -log10(p.adj))) +
  geom_point(aes(color = significance), alpha = 0.6, size = 2) +
  geom_text_repel(
    aes(label = label),
    max.overlaps = 15,
    size = 4,
    box.padding = 0.5
  ) +
  scale_color_manual(values = c("Significant" = "#E41A1C", "Not Significant" = "grey")) +
  labs(
    title = paste0(gOE, " In Silico Perturbation Analysis"),
    subtitle = "Overexpression Folds",
    x = "Z-score",
    y = "-log10(padj)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkgrey")
print(diff_reg_plot)
write.csv(dr_results, file = "RPS26_overexpression_results.csv")