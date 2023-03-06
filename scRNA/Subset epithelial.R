library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(dittoSeq)
library(pylr)
library(UCell)
library(Matrix)
library(ggpubr)
set.seed(123)

whole = readRDS('./whole.rds')
## subset Epithelium
epi = subset(whole, idents = c('AT2','Transitional','AT1','Secretory & Ciliated'))
DefaultAssay(epi) = 'RNA'
epi <- SCTransform(epi, verbose = FALSE)
epi <- RunPCA(epi, verbose = FALSE)

pct <- epi[["pca"]]@stdev / sum(epi[["pca"]]@stdev) * 100

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

epi <- RunUMAP(epi, dims = 1:19, verbose = FALSE)
epi <- FindNeighbors(epi, dims = 1:19, verbose = FALSE)
epi <- FindClusters(epi, verbose = FALSE, resolution = 0.5)
DimPlot(epi, reduction = "umap", label = T, label.size = 6, pt.size = 1.0)

## subset only Epcam+ cells
FeaturePlot(epi, reduction = "umap", features = c('Pecam1','Epcam','Col1a1','Ptprc','Msln'), order = TRUE, label = F, repel = T, label.size = 3, cols = c('gray','red'))
epi = subset(epi, idents = c('18','12','17','10','16','13'), invert = T)

epi <- SCTransform(epi, verbose = FALSE)
epi <- RunPCA(epi, verbose = FALSE)

pct <- epi[["pca"]]@stdev / sum(epi[["pca"]]@stdev) * 100

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

epi <- RunUMAP(epi, dims = 1:13, verbose = FALSE)
epi <- FindNeighbors(epi, dims = 1:13, verbose = FALSE)
epi <- FindClusters(epi, verbose = FALSE, resolution = 2.1)
DimPlot(epi, reduction = "umap", label = T, label.size = 6, pt.size = 1.0)

FeaturePlot(epi, reduction = "umap", features = c('Pecam1','Epcam','Col1a1','Ptprc','Msln'), order = TRUE, label = F, repel = T, label.size = 3, cols = c('gray','red'))

epi = subset(epi, idents = c('25'), invert = T)

epi <- SCTransform(epi, verbose = FALSE)
epi <- RunPCA(epi, verbose = FALSE)

pct <- epi[["pca"]]@stdev / sum(epi[["pca"]]@stdev) * 100

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

epi <- RunUMAP(epi, dims = 1:14, verbose = FALSE)
epi <- FindNeighbors(epi, dims = 1:14, verbose = FALSE)
epi <- FindClusters(epi, verbose = FALSE, resolution = 0.2)
epi <- FindSubCluster(epi, cluster = '5', graph.name='SCT_snn', subcluster.name='epi.cluster', resolution = 0.1)
Idents(epi) = epi@meta.data$epi.cluster

##Celltype annotation
epi.marker <- FindAllMarkers(epi, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
epi@meta.data$epi_celltype = Idents(epi) %>% as.character
epi@meta.data$epi_celltype = factor(epi@meta.data$epi_celltype, levels = c('AT1','AT2','AT2-Lyz1+','Transitional','Transitional-Lyz1+','Secretory','CiliaSecretory','Ciliated'))

epi = RenameIdents(epi,
'0' = 'AT2', 
'1' = 'Transitional',
'2' = 'AT2',
'3' = 'AT2-Lyz1+',
'4' = 'AT1',
'5_0' = 'Ciliated',
'5_1' = 'Secretory',
'5_2' = 'CiliaSecretory',
'6' = 'AT2',
'7' = 'Transitional-Lyz1+'
)

saveRDS(epi, '../epi.rds')


##Epithelium
p = DimPlot(epi, cols=c('#008A69','#FBBD80','#E59560','#BDD7EE','#5DA6E8', '#F58F9D','#3D49C8', '#A04AA4', '#808080'), group = 'epi_celltype')
LabelClusters(p,family = 'Arial', repel = T, size = 5, id = 'epi_celltype') +
theme(plot.title = element_text(family = 'Arial'), legend.text = element_text(size=20,family = 'Arial')) + 
theme(axis.title.x = element_text(family = 'Arial'), axis.title.y = element_text(family = 'Arial')) + 
theme(axis.text.x = element_text( family = 'Arial'), axis.text.y = element_text(family = 'Arial')) +
theme(legend.title = element_blank()) + guides(color = guide_legend(nrow = 4, byrow = T, override.aes = list(size=12))) +
theme(legend.position = 'bottom')

##Epithelium composition chart
ggplot(epi@meta.data, aes(x=state, fill=epi_celltype)) + geom_bar(colour = 'black', position = 'fill') + theme(axis.text.x = element_text(angle = 60, hjust = 1, size=20,color="black",family = 'Arial')) + 
scale_fill_manual(values = c('#008A69','#FBBD80','#E59560','#BDD7EE','#5DA6E8', '#F58F9D','#3D49C8', '#A04AA4'), name = '') + 
ggtitle('Composition stack chart') + ylab('Proportion') + 
theme(plot.title = element_text(size=15,family = 'Arial'), legend.text = element_text(size=20,family = 'Arial')) + 
theme(axis.title.x = element_text(size=20,family = 'Arial'), axis.title.y = element_text(size=20,family = 'Arial')) + 
theme(axis.text.y = element_text(size=20,family = 'Arial')) +
theme(legend.title = element_blank()) + guides(colour = guide_legend(override.aes = list(size=4)))

##Epithelium dot plot
dittoDotPlot(epi, vars = c('Mki67','Cdkn1a','Krt8','Lyz1','Sftpc','Sftpa1','Aqp5','Hopx','Scgb3a2','Scgb1a1','Dynlrb2','Foxj1'), group='epi_celltype', size = 12) +
theme(plot.title = element_text(size=15,family = 'Arial'), legend.text = element_text(size=18,family = 'Arial'), legend.title = element_text(size=20, family = 'Arial') ) + 
theme(axis.title.x = element_text(size=20,family = 'Arial'), axis.title.y = element_text(size=20,family = 'Arial')) + 
theme(axis.text.y = element_text(size=20,family = 'Arial'), axis.text.x = element_text(size=20,family = 'Arial',face = 'italic'))

##Ndfus2 expression in Epithelium
VlnPlot(epi, features=c('Ndufs2'), split.by = 'state', split.plot = F, group.by='epi_celltype', cols=c('#CF7C89', '#94D8F6'), slot = 'data', assay = 'SCT') + 
theme(axis.text.x = element_text(angle = 40, hjust = 1, size=20,color="black", family = 'Arial'), axis.text.y = element_text(size = 20, family = 'Arial'), axis.title.y = element_text(size=20,family = 'Arial')) + 
theme(plot.title = element_text(size = 24, face = 'italic', family = 'Arial')) + theme(strip.text.x = element_text(size = 20, family = 'Arial'), legend.text = element_text(size=20,family = 'Arial')) +
xlab('')+ylab('Expression level')+ggtitle('Ndufs2')


##krt expression
p1 = FeaturePlot(epi, reduction = "umap", features = c('Krt8','Krt18','Krt7','Krt19'),  order = TRUE, label = F, combine = FALSE)
fix.sc <- scale_color_gradientn( colours = c('lightgrey', 'blue'),  limits = c(0, 5))
p2 <- lapply(p1, function (x) x + fix.sc +  theme(plot.title = element_text(face = 'italic',family = 'Arial')))
CombinePlots(p2)

##cell cycle
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
#s.genes convert to mouse gene (m.s.genes)
#g2m.genes convert to mouse gene (m.g2m.genes)
epi <- CellCycleScoring(epi, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = TRUE)

##Atf expression
p1 = FeaturePlot(epi, reduction = "umap", features = c('Atf3','Atf4','Atf5','Atf6'),  order = TRUE, label = F, combine = FALSE)
fix.sc <- scale_color_gradientn( colours = c('lightgrey', 'blue'),  limits = c(0, 3))
p2 <- lapply(p1, function (x) x + fix.sc +  theme(plot.title = element_text(face = 'italic',family = 'TT Arial'), legend.text=element_text(family = 'TT Arial'), axis.text.x = element_text(color="black",family = 'TT Arial'), axis.text.y = element_text(color="black",family = 'TT Arial')))
CombinePlots(p2)

##ISR enrichment score
signatures = list(ISR = isr.genes)
my.matrix = GetAssayData(object = epi, slot = "data", assay = 'SCT')
u.scores = ScoreSignatures_UCell(my.matrix, features = signatures)
epi = AddMetaData(epi, metadata = as.data.frame(u.scores))

FeaturePlot(epi, reduction = "umap", features = 'ISR_UCell',  order = TRUE) +
ggtitle('ISR enrichment score') + theme(plot.title=element_text(family = 'Arial', size = 20), legend.text=element_text(family = 'Arial', size = 15), axis.text.x = element_text(color="black",family = 'Arial'), axis.text.y = element_text(color="black",family = 'Arial'))

ggplot(epi@meta.data, aes(x = epi_celltype, y = ISR_UCell)) + geom_violin(aes(fill = epi_celltype), scale = "width") + 
geom_boxplot(width = 0.1) + theme_bw() + stat_compare_means(label.y = 0.32,label.x = 1, family = 'Arial',size = 5) +
theme(plot.title = element_text(size = 20,family = 'Arial'), legend.text = element_text(size = 20,family = 'Arial')) + 
theme(axis.title.x = element_text(size = 20,family = 'Arial'), axis.title.y = element_text(size = 20,family = 'Arial')) + 
theme(axis.text.x = element_text(angle = 40, hjust = 1,size = 20, family = 'Arial'), axis.text.y = element_text(size = 20,family = 'Arial')) + NoLegend() + ylab('ISR score')

#scRNA ISR genes heatmap
heatmap_epi = epi
epi.genes <- rownames(heatmap_epi)
heatmap_epi <- ScaleData(heatmap_epi, features = epi.genes)
DoHeatmap(heatmap_epi, features=isr.genes, size = 5, slot = 'scale.data', assay='SCT') +
theme(plot.title = element_text(size = 20,family = 'Arial'), legend.text = element_text(size = 20,family = 'Arial')) + 
theme(axis.title.x = element_text(size = 20,family = 'Arial'), axis.title.y = element_text(size = 20,family = 'Arial')) + 
theme(axis.text.x = element_text(size = 20, family = 'Arial'), axis.text.y = element_text(size = 20,family = 'Arial')) +
guides(color = guide_legend(override.aes = list(size=4)))




