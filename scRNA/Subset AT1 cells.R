library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

epi = readRDS('../epi.rds')
## subset AT1
at1 = subset(epi, idents = c('AT1'))
DefaultAssay(at1) = 'RNA'
at1 <- SCTransform(at1, verbose = FALSE)
at1 <- RunPCA(at1, verbose = FALSE)

pct <- at1[["pca"]]@stdev / sum(at1[["pca"]]@stdev) * 100

cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
pcs <- min(co1, co2)
plot_df <- data.frame(pct = pct, cumu = cumu, rank = 1:length(pct))
 
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) +
  geom_text() +
  geom_vline(xintercept = 90, color = "grey") +
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
  
at1 <- RunUMAP(at1, dims = 1:14, verbose = FALSE)
at1 <- FindNeighbors(at1, dims = 1:14, verbose = FALSE)
at1 <- FindClusters(at1, verbose = FALSE, resolution = 0.1)
DimPlot(at1, reduction = "umap", label = T, label.size = 6, pt.size = 1.0)

#at1 subclustering markers
at1.markers <- FindAllMarkers(at1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
at1.markers_top50 = at1.markers %>% group_by(cluster) %>% slice_max(n=50, order_by=avg_log2FC)

