## whole cluster analysis
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

SHH1 = Read10X('./SHH1/outs/filtered_feature_bc_matrix')
SHH2 = Read10X('./SHH2/outs/filtered_feature_bc_matrix')
SHH3 = Read10X('./SHH3/outs/filtered_feature_bc_matrix')
SHH4 = Read10X('./SHH4/outs/filtered_feature_bc_matrix')
SHH5 = Read10X('./SHH5/outs/filtered_feature_bc_matrix')
SHH6 = Read10X('./SHH6/outs/filtered_feature_bc_matrix')
SHH7 = Read10X('./SHH7/outs/filtered_feature_bc_matrix')
SHH8 = Read10X('./SHH8/outs/filtered_feature_bc_matrix')

SHH1 = CreateSeuratObject(counts = SHH1, project = "shh1", min.features = 200, min.cells = 3)
SHH2 = CreateSeuratObject(counts = SHH2, project = "shh2", min.features = 200, min.cells = 3)
SHH3 = CreateSeuratObject(counts = SHH3, project = "shh3", min.features = 200, min.cells = 3)
SHH4 = CreateSeuratObject(counts = SHH4, project = "shh4", min.features = 200, min.cells = 3)
SHH5 = CreateSeuratObject(counts = SHH5, project = "shh5", min.features = 200, min.cells = 3)
SHH6 = CreateSeuratObject(counts = SHH6, project = "shh6", min.features = 200, min.cells = 3)
SHH7 = CreateSeuratObject(counts = SHH7, project = "shh7", min.features = 200, min.cells = 3)
SHH8 = CreateSeuratObject(counts = SHH8, project = "shh8", min.features = 200, min.cells = 3)

##scrublet results into seurat obj(SHH1 - SHH8)

#SHH1
doublet.file <- "../SHH1/scrublet.doublet_scores_obs.txt"
doublet.df <- read.table(doublet.file , header = F, sep = "\t")
colnames(doublet.df) <- c("barcodes", "scores", "threshold")
doublet.df$Doublet <- "singlet"
doublet.df$Doublet[doublet.df$scores > doublet.df$threshold] <- "doublet"
doublet.df$barcodes <- do.call(rbind, strsplit(x = as.character(doublet.df$barcodes), split = "\\-"))
metadata = data.frame(barcodes = rownames(SHH1@meta.data), SHH1@meta.data)
metadata$barcodes = do.call(rbind, strsplit(rownames(metadata), split = '\\-'))
test <- join(metadata, doublet.df, by = "barcodes")
SHH1@meta.data$Doublet <- test$Doublet

#SHH2
doublet.file <- "../SHH2/scrublet.doublet_scores_obs.txt"
doublet.df <- read.table(doublet.file , header = F, sep = "\t")
colnames(doublet.df) <- c("barcodes", "scores", "threshold")
doublet.df$Doublet <- "singlet"
doublet.df$Doublet[doublet.df$scores > doublet.df$threshold] <- "doublet"
doublet.df$barcodes <- do.call(rbind, strsplit(x = as.character(doublet.df$barcodes), split = "\\-"))
metadata = data.frame(barcodes = rownames(SHH2@meta.data), SHH2@meta.data)
metadata$barcodes = do.call(rbind, strsplit(rownames(metadata), split = '\\-'))
test <- join(metadata, doublet.df, by = "barcodes")
SHH2@meta.data$Doublet <- test$Doublet

#SHH3
doublet.file <- "../SHH3/scrublet.doublet_scores_obs.txt"
doublet.df <- read.table(doublet.file , header = F, sep = "\t")
colnames(doublet.df) <- c("barcodes", "scores", "threshold")
doublet.df$Doublet <- "singlet"
doublet.df$Doublet[doublet.df$scores > doublet.df$threshold] <- "doublet"
doublet.df$barcodes <- do.call(rbind, strsplit(x = as.character(doublet.df$barcodes), split = "\\-"))
metadata = data.frame(barcodes = rownames(SHH3@meta.data), SHH3@meta.data)
metadata$barcodes = do.call(rbind, strsplit(rownames(metadata), split = '\\-'))
test <- join(metadata, doublet.df, by = "barcodes")
SHH3@meta.data$Doublet <- test$Doublet

#SHH4
doublet.file <- "../SHH4/scrublet.doublet_scores_obs.txt"
doublet.df <- read.table(doublet.file , header = F, sep = "\t")
colnames(doublet.df) <- c("barcodes", "scores", "threshold")
doublet.df$Doublet <- "singlet"
doublet.df$Doublet[doublet.df$scores > doublet.df$threshold] <- "doublet"
doublet.df$barcodes <- do.call(rbind, strsplit(x = as.character(doublet.df$barcodes), split = "\\-"))
metadata = data.frame(barcodes = rownames(SHH4@meta.data), SHH4@meta.data)
metadata$barcodes = do.call(rbind, strsplit(rownames(metadata), split = '\\-'))
test <- join(metadata, doublet.df, by = "barcodes")
SHH4@meta.data$Doublet <- test$Doublet

#SHH5
doublet.file <- "../SHH5/scrublet.doublet_scores_obs.txt"
doublet.df <- read.table(doublet.file , header = F, sep = "\t")
colnames(doublet.df) <- c("barcodes", "scores", "threshold")
doublet.df$Doublet <- "singlet"
doublet.df$Doublet[doublet.df$scores > doublet.df$threshold] <- "doublet"
doublet.df$barcodes <- do.call(rbind, strsplit(x = as.character(doublet.df$barcodes), split = "\\-"))
metadata = data.frame(barcodes = rownames(SHH5@meta.data), SHH5@meta.data)
metadata$barcodes = do.call(rbind, strsplit(rownames(metadata), split = '\\-'))
test <- join(metadata, doublet.df, by = "barcodes")
SHH5@meta.data$Doublet <- test$Doublet

#SHH6
doublet.file <- "../SHH6/scrublet.doublet_scores_obs.txt"
doublet.df <- read.table(doublet.file , header = F, sep = "\t")
colnames(doublet.df) <- c("barcodes", "scores", "threshold")
doublet.df$Doublet <- "singlet"
doublet.df$Doublet[doublet.df$scores > doublet.df$threshold] <- "doublet"
doublet.df$barcodes <- do.call(rbind, strsplit(x = as.character(doublet.df$barcodes), split = "\\-"))
metadata = data.frame(barcodes = rownames(SHH6@meta.data), SHH6@meta.data)
metadata$barcodes = do.call(rbind, strsplit(rownames(metadata), split = '\\-'))
test <- join(metadata, doublet.df, by = "barcodes")
SHH6@meta.data$Doublet <- test$Doublet

#SHH7
doublet.file <- "../SHH7/scrublet.doublet_scores_obs.txt"
doublet.df <- read.table(doublet.file , header = F, sep = "\t")
colnames(doublet.df) <- c("barcodes", "scores", "threshold")
doublet.df$Doublet <- "singlet"
doublet.df$Doublet[doublet.df$scores > doublet.df$threshold] <- "doublet"
doublet.df$barcodes <- do.call(rbind, strsplit(x = as.character(doublet.df$barcodes), split = "\\-"))
metadata = data.frame(barcodes = rownames(SHH7@meta.data), SHH7@meta.data)
metadata$barcodes = do.call(rbind, strsplit(rownames(metadata), split = '\\-'))
test <- join(metadata, doublet.df, by = "barcodes")
SHH7@meta.data$Doublet <- test$Doublet

#SHH8
doublet.file <- "../SHH8/scrublet.doublet_scores_obs.txt"
doublet.df <- read.table(doublet.file , header = F, sep = "\t")
colnames(doublet.df) <- c("barcodes", "scores", "threshold")
doublet.df$Doublet <- "singlet"
doublet.df$Doublet[doublet.df$scores > doublet.df$threshold] <- "doublet"
doublet.df$barcodes <- do.call(rbind, strsplit(x = as.character(doublet.df$barcodes), split = "\\-"))
metadata = data.frame(barcodes = rownames(SHH8@meta.data), SHH8@meta.data)
metadata$barcodes = do.call(rbind, strsplit(rownames(metadata), split = '\\-'))
test <- join(metadata, doublet.df, by = "barcodes")
SHH8@meta.data$Doublet <- test$Doublet

merge12 = merge(SHH1, SHH2, add.cell.id = c('shh1','shh2'))
merge34 = merge(SHH3, SHH4, add.cell.id = c('shh3','shh4'))
merge56 = merge(SHH5, SHH6, add.cell.id = c('shh5','shh6'))
merge78 = merge(SHH7, SHH8, add.cell.id = c('shh7','shh8'))
merge1234 = merge(merge12, merge34)
merge5678 = merge(merge56, merge78)
merge = merge(merge1234, merge5678)

merge[["percent.mt"]] <- PercentageFeatureSet(merge, pattern = "^mt-")

merge@meta.data$state = merge@meta.data$orig.ident %>% as.character
merge@meta.data$state[merge@meta.data$state == 'SHH1'] = 'NDUFS2 cKO'
merge@meta.data$state[merge@meta.data$state == 'SHH2'] = 'NDUFS2 cKO'
merge@meta.data$state[merge@meta.data$state == 'SHH5'] = 'NDUFS2 cKO'
merge@meta.data$state[merge@meta.data$state == 'SHH8'] = 'NDUFS2 cKO'
merge@meta.data$state[merge@meta.data$state == 'SHH3'] = 'NDUFS2 control'
merge@meta.data$state[merge@meta.data$state == 'SHH4'] = 'NDUFS2 control'
merge@meta.data$state[merge@meta.data$state == 'SHH6'] = 'NDUFS2 control'
merge@meta.data$state[merge@meta.data$state == 'SHH7'] = 'NDUFS2 control'

merge@meta.data$sample = merge@meta.data$orig.ident %>% as.character
merge@meta.data$sample[merge@meta.data$sample == 'SHH1'] = 'NDUFS2 cKO M1'
merge@meta.data$sample[merge@meta.data$sample == 'SHH2'] = 'NDUFS2 cKO F1'
merge@meta.data$sample[merge@meta.data$sample == 'SHH5'] = 'NDUFS2 cKO F2'
merge@meta.data$sample[merge@meta.data$sample == 'SHH8'] = 'NDUFS2 cKO M2'
merge@meta.data$sample[merge@meta.data$sample == 'SHH3'] = 'NDUFS2 control F1'
merge@meta.data$sample[merge@meta.data$sample == 'SHH4'] = 'NDUFS2 control F2'
merge@meta.data$sample[merge@meta.data$sample == 'SHH6'] = 'NDUFS2 control M1'
merge@meta.data$sample[merge@meta.data$sample == 'SHH7'] = 'NDUFS2 control M2'

##Doublet 2409, Singlet 60657
merge <- subset(x = merge, subset = Doublet == "singlet")

##visualize nFeatrue_RNA & nCount_RNA & percent.mt
VlnPlot(merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

merge_filter = subset(merge, subset = nFeature_RNA > 500 & percent.mt < 25)

whole <- SCTransform(merge_filter, verbose = FALSE)
whole <- RunPCA(whole, verbose = FALSE)

## dims
pct <- whole[["pca"]]@stdev / sum(whole[["pca"]]@stdev) * 100

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

whole <- RunUMAP(whole, dims = 1:22, verbose = FALSE)
whole <- FindNeighbors(whole, dims = 1:22, verbose = FALSE)
whole <- FindClusters(whole, verbose = FALSE, resolution = 0.6)
DimPlot(whole, reduction = "umap", label = T, label.size = 6, repel = T)

##Celltype annotation
whole.marker <- FindAllMarkers(whole, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

whole = RenameIdents(whole,
'0' = 'AT2', 
'1' = 'B',
'2' = 'gCap',
'3' = 'T',
'4' = 'B',
'5' = 'FIB.1',
'6' = 'MP',
'7' = 'Transitional',
'8' = 'gCap',
'9' = 'Neutrophils',
'10' = 'Art',
'11' = 'CM',
'12' = 'NCM',
'13' = 'aCap',
'14' = 'IM',
'15' = 'DC.2',
'16' = 'SM.2',
'17' ='FIB.2',
'18' ='ILC',
'19' ='AT1',
'20' ='SM.1',
'21' ='NK',
'22' ='CD8+T',
'23' ='DC.1',
'24'='Secretory & Ciliated',
'25'='Neutrophils',
'26'='Treg',
'27'='Pl.DC',
'28'='Mesothelium',
'29'='Lymph EC',
'30'='Pl.MP',
'31'='Ccr7+DC'
)

whole@meta.data$chart = Idents(whole) %>% as.character

saveRDS(whole, '../whole.rds')

## whole clusters
p = DimPlot(whole)
LabelClusters(p,family = 'Arial', repel = T, size = 5, id = 'ident') +
theme(plot.title = element_text(family = 'Arial'), legend.text = element_text(family = 'Arial')) + 
theme(axis.title.x = element_text(family = 'Arial'), axis.title.y = element_text(family = 'Arial')) + 
theme(axis.text.x = element_text( family = 'Arial'), axis.text.y = element_text(family = 'Arial')) +
theme(legend.title = element_blank()) + NoLegend()

## whole clusters split.by state
DimPlot(whole, group = 'state', cols=c('#CF7C89', '#94D8F6')) +
theme(plot.title = element_text(family = 'Arial'), legend.text = element_text(size=20, family = 'Arial')) + 
theme(axis.title.x = element_text(family = 'Arial'), axis.title.y = element_text(family = 'Arial')) + 
theme(axis.text.x = element_text( family = 'Arial'), axis.text.y = element_text(family = 'Arial')) +
theme(legend.title = element_blank()) + theme(legend.position = c(0.02,0.95)) + guides(colour = guide_legend(override.aes = list(size=4))) + ggtitle('')

## tdTomato expression
FeaturePlot(whole, reduction = "umap", features = c('tdTomato'),  order = TRUE, label = F) +  theme(plot.title = element_text(face = 'italic',family = 'Arial'))

## Ndfus2 expression in whole
VlnPlot(whole, features=c('Ndufs2'), split.by = 'state', split.plot = F, group.by='chart', cols=c('#CF7C89', '#94D8F6'), slot = 'data', assay = 'SCT') + 
theme(axis.text.x = element_text(angle = 40, hjust = 1, size=20,color="black", family = 'Arial'), axis.text.y = element_text(size = 20, family = 'Arial'), axis.title.y = element_text(size=20,family = 'Arial')) + 
theme(plot.title = element_text(size = 24, face = 'italic', family = 'Arial')) + theme(strip.text.x = element_text(size = 20, family = 'Arial'), legend.text = element_text(size=20,family = 'Arial')) +
xlab('')+ylab('Expression level')+ggtitle('Ndufs2')

## whole cluster ISR
signatures = list(ISR = isr.genes)
my.matrix = GetAssayData(object = whole, slot = "data", assay = 'SCT')
u.scores = ScoreSignatures_UCell(my.matrix, features = signatures)
whole = AddMetaData(whole, metadata = as.data.frame(u.scores))

##ISR enrichment score
ggplot(whole@meta.data, aes(x = chart, y = ISR_UCell)) + geom_violin(aes(fill = chart), scale = "width") + 
geom_boxplot(width = 0.1) + theme_bw() + stat_compare_means(label.y = 0.32,label.x = 1, family = 'Arial',size = 5) +
theme(plot.title = element_text(size = 20,family = 'Arial'), legend.text = element_text(size = 20,family = 'Arial')) + 
theme(axis.title.x = element_text(size = 20,family = 'Arial'), axis.title.y = element_text(size = 20,family = 'Arial')) + 
theme(axis.text.x = element_text(angle = 40, hjust = 1,size = 20, family = 'Arial'), axis.text.y = element_text(size = 20,family = 'Arial')) + NoLegend() + ylab('ISR score')

## subset Mki67+ cells
my.matrix = GetAssayData(object = whole, slot = "data", assay='SCT')
my.matrix = my.matrix %>% t %>% as.data.frame
mki_cells = filter(my.matrix, Mki67 > 0) %>% rownames()


## whole cluster heatmap
whole.genes <- rownames(whole)
whole <- ScaleData(whole, features = whole.genes)
marker = c('Sftpc','Sftpa1','Sftpb','Lamp3','Sftpd','Lyz2',
           'Krt8','Krt18','Atf5','Cdkn1a',
           'Hopx','Aqp5','Vegfa','Col4a3',
           'Scgb1a1','Scgb3a2','Dynlrb2','Foxj1','Tubb4b',
           'Inmt','Col3a1','Pdgfra',
           'Col1a1','Col1a2','Timp1','Hhip','Sfrp1',
           'Gucy1a1','Cox4i2','Pdgfrb',
           'Eln','Acta2','Tagln','Myh11',
           'Upk3b','Sfrp2','Msln',
           'Cd93','Ptprb','Gpihbp1','Kit',
           'Ednrb','Kdr','Car4',
           'Cxcl12','Vwf',
           'Ccl21a','Mmrn1','Nrp2','Flt4','Prss23',
           'Igkc','Cd79a','Ly6d','Ms4a1','Cd74',
           'Ms4a4b','Cd3d','Cd3g',
           'Cd8b1','Ly6c2',
           'Trbc2','Trbc1','S100a10','AW112010',
           'S100a9','S100a8','Retnlg','Ifitm1','Csf3r','Cxcr2',
           'Plac8','Ifitm3','S100a4','Ifi27l2a','Ccr2',
           'Pou2f2','Cx3cr1','Csf1r',
           'Chil3','Ctsd','Plet1','Lpl','Atp6v0d2',
           'C1qa','C1qb','C1qc','Apoe','Pf4',
           'Mki67','Top2a',
           'Cst3','H2-Ab1','H2-Eb1','H2-Aa',
           'Mgl2','Cd209a',
           'Ccl5','Fscn1','Ccr7','Traf1',
           'Gzma','Nkg7','Klrd1',
           'Il7r','Cxcr6')
color = rep(c('#DAAADB','#C785C8','#7A297B','#A349A4','#7030A0','#511252','#0010FF','#8389E0','#3F48CC','#232B99','#101566','#60C5F1','#00A2E8','#0070C0','#007AAE','#005174','#6BD089','#52AC6D','#138535','#085820','#FFB27D','#9F9700','#BFB500','#A89F00','#807900','#F8A1A4','#F47378','#B21016','#77070B'),c(2,3,4,0,2,4,2,5,5,3,5,6,4,2,3,5,5,2,3,4,3,4,3,5,3,5,4,4,6))
DoHeatmap(whole, features=marker, slot = 'scale.data', assay='SCT', group.by='chart', group.colors=c('#77070B','#B21016','#F47378','#F8A1A4','#807900','#A89F00','#BFB500','#9F9700','#FFB27D','#085820','#138535','#52AC6D','#6BD089','#005174','#007AAE','#0070C0','#00A2E8','#60C5F1','#101566','#232B99','#3F48CC','#8389E0','#0010FF','#511252','#7030A0','#A349A4','#7A297B','#C785C8','#DAAADB')) +
theme(plot.title = element_text(size = 20,family = 'Arial'), legend.text = element_text(size = 20,family = 'Arial')) + 
theme(axis.title.y = element_text(family = 'Arial')) + 
theme(axis.text.y = element_text(family = 'Arial',color=color))
ggsave(path = './', device='tiff', dpi=1000, filename='heatmap.tiff')
dev.off()






