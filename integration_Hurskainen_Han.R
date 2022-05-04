## scanpy UMAP
sample.info = read.csv('../integ_Hur_Han/sample.info.csv')
umap = read.csv('../integ_Hur_Han/umap.csv')
## we have to match the cell ids of seurat obj. and scanpy obj. ##

umap = umap %>% as.data.frame()
umap = umap[,-1]
colnames(umap) = c('UMAP_1','UMAP_2')
rowanems(umap) = sample.info$X


## load Hurskainen et al. epithelium seurat obj
hur_epi = readRDS('../hur_epi.rds')
hur_epi@meta.data$projcet = 'Hurskanien et al.'
hur_epi@meta.data$epi_celltype = Idents(hur_epi)

## load epithelium seurat obj from current work
epi = readRDS('../epi.rds')
epi@meta.data$Age = 'P21'
epi@meta.data$project = 'Han et al.'
epi@meta.data$epi_celltype = Idents(epi)

## merge data
seurat.list <- SplitObject(all, split.by="project")
for (i in 1:length(seurat.list)) {
    seurat.list[[i]] <- SCTransform(seurat.list[[i]], 
                                    vars.to.regress=c("percent.mt"),
                                    verbose = FALSE)
}
seurat.features <- SelectIntegrationFeatures(seurat.list, nfeatures = 3000)
seurat.list <- PrepSCTIntegration(seurat.list, anchor.features = seurat.features, verbose = T)
seurat.anchors <- FindIntegrationAnchors(seurat.list, normalization.method = "SCT", anchor.features = seurat.features, verbose = T)
seurat <- IntegrateData(seurat.anchors, normalization.method = "SCT",  verbose = T)

seurat <- RunPCA(seurat, verbose=F)
seurat <- RunUMAP(seurat, dims=1:40)

Idents(seurat) = seurat@meta.data$epi_celltype
seurat@meta.data$state[seurat@meta.data$Oxygen == 'Hyperoxia'] = 'Hyperoxia'
seurat@meta.data$state[seurat@meta.data$Oxygen == 'Normoxia'] = 'Normoxia'

## extract info to use in scanpy
# Cell ids
write.csv(Cells(seurat), file = "../integ_Hur_Han/cellID_obs.csv", row.names = FALSE)
# Clusters
write.csv(select(seurat@meta.data, epi_celltype), '../integ_Hur_Han/clusters.csv')
# Cluster colors
ident.colors <- (scales::hue_pal())(n = length(x = levels(x = seurat)))
names(x = ident.colors) <- levels(x = seurat)
cell.colors <- ident.colors[Idents(object = seurat)]
names(x = cell.colors) <- colnames(x = seurat)
write.csv(cell.colors, '../integ_Hur_Han/cell.colors.csv')
# Age
write.csv(select(seurat@meta.data, Age), '../integ_Hur_Han/Age.csv')
# project
write.csv(select(seurat@meta.data, project), '../integ_Hur_Han/project.csv')


## UMAP embedding
seurat[['umap']]@cell.embeddings = umap %>% as.matrix

## Dim plot
p = DimPlot(seurat, group = 'epi_celltype')
LabelClusters(p,family = 'Arial', repel = T, size = 5, id = 'epi_celltype') +
theme(plot.title = element_text(family = 'Arial'), legend.text = element_text(size = 30,family = 'Arial')) + 
theme(axis.title.x = element_text(family = 'Arial'), axis.title.y = element_text(family = 'Arial')) + 
theme(axis.text.x = element_text( family = 'Arial'), axis.text.y = element_text(family = 'Arial')) +
theme(legend.title = element_blank()) + ggtitle('') + theme(legend.position = 'bottom') + guides(color = guide_legend(nrow = 4, byrow = T, override.aes = list(size=12)))


## calculate ISR score
signatures = list(ISR = isr.genes)
my.matrix = GetAssayData(object = seurat, slot = "data", assay = 'SCT')
u.scores = ScoreSignatures_UCell(my.matrix, features = signatures)
seurat = AddMetaData(seurat, metadata = as.data.frame(u.scores))



#subset NDUFS2 cKO & Hyperoxia
seurat@meta.data$tmp = paste0(seurat@meta.data$state,'_',seurat@meta.data$epi_celltype)
seurat_tran_hyat2 = subset(seurat, idents = c('NDUFS2 cKO_Transitional','Hyperoxia_AT2'))
ggplot(seurat_tran_hyat2@meta.data, aes(x = project, y = ISR_UCell)) + 
geom_violin(aes(fill = project), scale = "width") + geom_boxplot(width = 0.1) + theme_bw() + stat_compare_means(label.y = 0.35, label.x = 0.6, size = 5, family = 'Arial') +
theme(plot.title = element_text(size = 20,family = 'Arial'), legend.text = element_text(size = 20,family = 'Arial')) + 
theme(axis.title.x = element_text(size = 20,family = 'Arial'), axis.title.y = element_text(size = 20,family = 'Arial')) + 
theme(axis.text.x = element_text(size = 20,family = 'Arial'), axis.text.y = element_text(size = 20,family = 'Arial')) +
theme(legend.title = element_blank()) + NoLegend() + ggtitle('') + ylab('ISR Score') + xlab('')




