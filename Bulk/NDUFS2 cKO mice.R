library(edgeR)
library(tidyverse)
library(ggrepel)
library(pheatmap)
library(viridis)
library(dplyr)
library(data.table)

rawdata=read.table("./ndufs2_bulk/gene_read_counts_table_all_final.tsv", header=TRUE, stringsAsFactors=FALSE, row.names=1)
sample_info = read.csv("./ndufs2_bulk/sample_info.csv")

sample_info$Group[sample_info$Group == 'con_nd'] = 'NDUFS2 control'
sample_info$Group[sample_info$Group == 'con_nd_ndi+'] = 'NDUFS2 control/NDI1'
sample_info$Group[sample_info$Group == 'ko_nd'] = 'NDUFS2 cKO'
sample_info$Group[sample_info$Group == 'ko_nd_ndi+'] = 'NDUFS2 cKO/NDI1'

class = factor(sample_info$Group, levels = c('NDUFS2 control','NDUFS2 control/NDI1','NDUFS2 cKO','NDUFS2 cKO/NDI1'))
sex = factor(sample_info$Sex, levels = c('F', 'M'))

genes=rownames(rawdata)

y <- DGEList(counts=rawdata, genes=genes, group=class)
y_filtered_norm <- calcNormFactors(y[filterByExpr(y),keep.lib.sizes = FALSE])

dge_design <- model.matrix(~sex + class)
y_dsp <- estimateDisp(y_filtered_norm, dge_design, robust = TRUE)
y_qlf_test <- glmQLFit(y_dsp, dge_design, robust = TRUE)
y_qlf_test <- glmQLFTest(y_qlf_test)

sample_info$Sex = factor(sample_info$Sex, levels = c('F','M'))
sample_info$Group = factor(sample_info$Group, levels = c('NDUFS2 control','NDUFS2 control/NDI1','NDUFS2 cKO','NDUFS2 cKO/NDI1'))

reg_sum <- decideTests(y_qlf_test)
group_annot <- sample_info %>% arrange(Group)
group_annot = group_annot[order(group_annot$Group, group_annot$Sex),]

rownames(group_annot) = group_annot$Sample.name
group_annot = group_annot[,-1]

ann_colors = list(Group = c('NDUFS2 control' = '#F8766D', 'NDUFS2 control/NDI1' = '#00BFC4','NDUFS2 cKO' = '#7CAE00', 'NDUFS2 cKO/NDI1' = '#C77CFF'), Sex = c('F' = '#818181', 'M' = '#C0BEBF'))

isr_sum <- decideTests(y_qlf_test, p.value = 1)

isr_counts <- as.logical(abs(isr_sum)) %>%
  y_filtered_norm[.,] %>%
  edgeR::cpm(.) %>%
  as.data.frame() %>% 
  dplyr::select(rownames(group_annot))
  
Atf=isr_counts[c('Atf3','Atf4','Atf5','Atf6','Ddit3'),]
newnames <- lapply(rownames(Atf),function(x) bquote(italic(.(x))))   
colors <- rev(colorRampPalette(RColorBrewer::brewer.pal(10, "RdBu"))(256))
pheatmap(isr_counts[c('Atf3','Atf4','Atf5','Atf6','Ddit3'),],
         annotation_col = dplyr::select(group_annot, c(Group,Sex)),
         annotation_colors = ann_colors,
         cluster_cols = F,
         cluster_rows = T,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         show_rownames = T,
         show_colnames = F,
         fontsize = 20,
         family = 'Arial',
         border_color = NA,
         #breaks = seq(from=-2,to=2,length.out=100),
         labels_row = as.expression(newnames),
         cellwidth = 20, cellheight = 40,
         gaps_col = c(8,15,22,30),
         scale = 'row',
         color = colors)


isr = read.table('.//isr.txt')
isr = isr$V1 %>% unique()
isr_heatmap = isr_counts[isr,]
isr_heatmap = na.omit(isr_heatmap)
newnames <- lapply(isr,function(x) bquote(italic(.(x))))  
pheatmap(isr_heatmap,
         annotation_col = dplyr::select(group_annot, c(Group,Sex)),
         annotation_colors = ann_colors,
         cluster_cols = F,
         cluster_rows = T,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         show_rownames = F,
         show_colnames = F,
         fontsize = 14,
         border_color = NA,
         fontsize_row = 15,
         family = 'Arial',
         #breaks = seq(from=-2,to=2,length.out=100),
         #labels_row = as.expression(newnames),
         #cellwidth = 20, cellheight = 40,
         gaps_col = c(8,15,22,30),
         scale = 'row',
         color = colors)
  
