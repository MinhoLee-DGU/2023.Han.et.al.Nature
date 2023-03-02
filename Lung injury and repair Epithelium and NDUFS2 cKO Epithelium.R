library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(dittoSeq)
library(UCell)
library(Matrix)
library(ggpubr)
set.seed(123)

## adult bleomycin injury model, GSE141259, Strunz et al.
load('./EpiHiRes_seurat.RData')
strunz = UpdateSeuratObject(subset)

strunz_meta = fread('./GSE141259_HighResolution_cellinfo.csv')
strunz_meta = strunz_meta %>% as.data.frame()
rownames(strunz_meta) = strunz_meta$cell_barcode

strunz_epi = subset(strunz, cells=rownames(strunz_meta))
strunz_epi = AddMetaData(strunz_epi, select(strunz_meta, cell_type))

strunz_umap = select(strunz_meta, c(umap_1,umap_2))
colnames(strunz_umap) = c('UMAP_1','UMAP_2')

strunz_epi[['umap']]@cell.embeddings = strunz_umap

han_epi = readRDS('./han_epi.rds')

DefaultAssay(han_epi) = 'RNA'
DefaultAssay(strunz_epi) = 'RNA'

strunz_epi[["percent.mt"]] <- PercentageFeatureSet(strunz_epi, pattern = "^mt-")
strunz_epi@meta.data$project = 'Strunz'
strunz_epi@meta.data$chart = strunz_epi@meta.data$cell_type

han_epi@meta.data$project = 'Han et al.'
strunz_epi@meta.data$project = 'Strunz et al.'

all = merge(strunz_epi, han_epi)
DefaultAssay(all) = 'RNA'
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
rm(seurat.anchors, seurat.list, seurat.features)

seurat <- RunPCA(seurat, verbose=F)
seurat <- RunUMAP(seurat, dims=1:14)

signatures = list(ISR = isr.genes)
my.matrix = GetAssayData(object = seurat, slot = "data", assay = 'SCT')
u.scores = ScoreSignatures_UCell(my.matrix, features = signatures)
seurat = AddMetaData(seurat, metadata = as.data.frame(u.scores))




## adult bleomycin injury model, GSE145031, Choi et al.
choi_epi = readRDS('./choi_epi.rds')
han_epi = readRDS('./han_epi.rds')

choi_epi@meta.data$project = 'Choi et al.'
han_epi@meta.data$project = 'Han et al.'

all = merge(choi_epi, han_epi)
DefaultAssay(all) = 'RNA'
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
rm(seurat.anchors, seurat.list, seurat.features)

seurat <- RunPCA(seurat, verbose=F)
seurat <- RunUMAP(seurat, dims=1:15)

signatures = list(ISR = isr.genes)
my.matrix = GetAssayData(object = seurat, slot = "data", assay = 'SCT')
u.scores = ScoreSignatures_UCell(my.matrix, features = signatures)
seurat = AddMetaData(seurat, metadata = as.data.frame(u.scores))


## mosue lung organoids, GSE141634, Kobayashi et al.
koba_epi = readRDS('./koba_epi.rds')
han_epi = readRDS('./han_epi.rds')

koba_epi@meta.data$project = 'Kobayashi et al.'
han_epi@meta.data$project = 'Han et al.'

all = merge(koba_epi, han_epi)
DefaultAssay(all) = 'RNA'
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
rm(seurat.anchors, seurat.list, seurat.features)

seurat <- RunPCA(seurat, verbose=F)
seurat <- RunUMAP(seurat, dims=1:15)

signatures = list(ISR = isr.genes)
my.matrix = GetAssayData(object = seurat, slot = "data", assay = 'SCT')
u.scores = ScoreSignatures_UCell(my.matrix, features = signatures)
seurat = AddMetaData(seurat, metadata = as.data.frame(u.scores))
