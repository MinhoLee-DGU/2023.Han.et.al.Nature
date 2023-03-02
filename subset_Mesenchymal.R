library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(dittoSeq)

whole = readRDS('./whole.rds')

mes = subset(whole, idents = c('SM.2','SM.1','FIB.1','FIB.2'))
mes <- SCTransform(mes, verbose = FALSE)
mes <- RunPCA(mes, verbose = FALSE)
mes <- RunUMAP(mes, dims = 1:20, verbose = FALSE)
mes <- FindNeighbors(mes, dims = 1:20, verbose = FALSE)
mes <- FindClusters(mes, verbose = FALSE, resolution = 0.2)
DimPlot(mes, reduction = "umap", label = T, label.size = 6, pt.size = 1.0)

FeaturePlot(mes, reduction = "umap", features = c('Epcam','Pecam1','Ptprc','Msln','Col1a1'), order = TRUE, label = F)
mes = subset(mes, idents= c('8','7'), invert = T)
DefaultAssay(mes) = 'RNA'

mes <- SCTransform(mes, verbose = FALSE)
mes <- RunPCA(mes, verbose = FALSE)
mes <- RunUMAP(mes, dims = 1:15, verbose = FALSE)
mes <- FindNeighbors(mes, dims = 1:15, verbose = FALSE)
mes <- FindClusters(mes, verbose = FALSE, resolution = 0.5)
DimPlot(mes, reduction = "umap", label = T, label.size = 6, pt.size = 1.0)


han_mes = RenameIdents(mes,
'0' = 'Distal', 
'1' = 'Pericytes',
'2' = 'Distal',
'3' = 'Distal',
'4' = 'Vascular SM',
'5' = 'Distal',
'6' = 'Distal-Sfrp1+',
'7' = 'Proximal',
'8' = 'Airway SM',
'9' = 'Myofib',
'10' = 'Myofib-Sfrp1+'
)


#heatmap
mes.genes <- rownames(han_mes)
han_mes <- ScaleData(han_mes, features = mes.genes)
mes_marker = c('Postn','Cox4i2','Gucy1a1','Pdgfrb','Notch3','Acta2','Tagln','Myh11','Actc1', 'Pi16','Serpinf1','Clec3b','Ccl11','Twist2','Wnt2','Tcf21','Npnt','Pdgfra','Vegfa','Bmp3','Gyg','Ces1d','Timp1','Sfrp1','Runx1','Chl1','Cp','Hhip','Cdh4','Wnt5a','Fgf18','Aspn','Lgr5','Lgr6','Tgfbi')
color = rep(c('#C785C8','#7A297B','#60C5F1','#007AAE','#F58F9D','#E59560','#FBBD80','#008A69'),c(1,7,5,8,5,1,3,5))
DoHeatmap(han_mes, features=mes_marker, slot = 'scale.data', assay='SCT', group.by='chart', group.colors=c('#008A69','#FBBD80','#E59560','#F58F9D','#007AAE','#60C5F1','#7A297B', '#C785C8')) +
theme(plot.title = element_text(size = 20,family = 'Arial'), legend.text = element_text(size = 15,family = 'Arial')) + 
theme(axis.title.y = element_text(family = 'Arial')) + 
theme(axis.text.y = element_text(family = 'Arial',color=color))

#UMAP
p = DimPlot(han_mes, cols=c('#008A69','#FBBD80','#E59560','#F58F9D','#007AAE','#60C5F1','#7A297B', '#C785C8'), group = 'chart')
LabelClusters(p,family = 'Arial', repel = T, size = 5, id = 'chart') +
theme(plot.title = element_text(family = 'Arial'), legend.text = element_text(size=20,family = 'Arial')) + 
theme(axis.title.x = element_text(family = 'Arial'), axis.title.y = element_text(family = 'Arial')) + 
theme(axis.text.x = element_text( family = 'Arial'), axis.text.y = element_text(family = 'Arial')) +
theme(legend.title = element_blank()) + guides(color = guide_legend(nrow = 4, byrow = T, family='Arial',override.aes = list(size=12))) +
theme(legend.position = 'bottom')

#composition chart
ggplot(han_mes@meta.data, aes(x=state, fill=chart)) + geom_bar(colour = 'black', position = 'fill') + theme(axis.text.x = element_text(angle = 60, hjust = 1, size=20,color="black",family = 'Arial')) + 
scale_fill_manual(values = c('#008A69','#FBBD80','#E59560','#F58F9D','#007AAE','#60C5F1','#7A297B', '#C785C8'), name = '') + 
ggtitle('Composition stack chart') + ylab('Proportion') + 
theme(plot.title = element_text(size=15,family = 'Arial'), legend.text = element_text(size=20,family = 'Arial')) + 
theme(axis.title.x = element_text(size=20,family = 'Arial'), axis.title.y = element_text(size=20,family = 'Arial')) + 
theme(axis.text.y = element_text(size=20,family = 'Arial')) +
theme(legend.title = element_blank()) + guides(colour = guide_legend(override.aes = list(size=4)))

#FeaturePlot
FeaturePlot(han_mes, reduction = "umap", features = c('Pdgfra','Pdgfrb','Acta2','Actc1','Wnt5a','Wnt2','Fgf18','Tgfbi','Hhip','Pi16','Timp1','Runx1','Sfrp1'),ncol=3, order = TRUE) + theme(plot.title=element_text(family = 'Arial', size = 20), legend.text=element_text(family = 'Arial', size = 15), axis.text.x = element_text(color="black",family = 'Arial'), axis.text.y = element_text(color="black",family = 'Arial'))










