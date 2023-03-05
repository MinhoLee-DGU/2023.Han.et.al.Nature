library(dplyr)
library(data.table)
library(tidyverse)
library(ggrepel)
library(edgeR)
library(pheatmap)
library(viridis)


dir = './p6_wt_bulk'
sample_info = read.csv('./p6_wt_bulk/sample_info_p6_wt.csv')
rawdata=read.table("./p6_wt_bulk/gene_read_counts_table_all_final.tsv", header=TRUE, stringsAsFactors=FALSE, row.names=1)

class = factor(sample_info$Group, levels = c('WT','P'))
sample_info$cul = sample_info$cul %>% as.character()
cul = factor(sample_info$cul, levels = c('0', '72'))

genes=rownames(rawdata)

y <- DGEList(counts=rawdata, genes=genes, group=class)
y_filtered_norm <- calcNormFactors(y[filterByExpr(y),keep.lib.sizes = FALSE])

dge_design <- model.matrix(~cul + class)
y_dsp <- estimateDisp(y_filtered_norm, dge_design, robust = TRUE)
y_qlf_test <- glmQLFit(y_dsp, dge_design, robust = TRUE)
y_qlf_test <- glmQLFTest(y_qlf_test)

sample_info$Sex = factor(sample_info$Sex, levels = c('F','M'))
sample_info$Group = factor(sample_info$Group, levels = c('WT','P'))
sample_info$ori = sample_info$ori %>% as.character()
sample_info$ori = factor(sample_info$ori, levels = c('2','3','1','4'))

reg_sum <- decideTests(y_qlf_test)
group_annot <- sample_info %>% arrange(Group)
group_annot = group_annot[order(group_annot$ori, group_annot$Group),]

rownames(group_annot) = group_annot$Sample.name
group_annot = group_annot[,-1]

group_annot$cultured = group_annot$cultured %>% as.character
colnames(group_annot)[4] = 'state'
ann_colors = list(Group = c('WT' = '#F8766D', 'P' = '#00BFC4'), state = c( '0' = '#C0BEBF', '72' = '#818181'))

isr_sum <- decideTests(y_qlf_test, p.value = 1)

isr_counts <- as.logical(abs(isr_sum)) %>%
  y_filtered_norm[.,] %>%
  edgeR::cpm(.) %>%
  as.data.frame() %>% 
  dplyr::select(rownames(group_annot))
  
scmar = isr_counts[c('Mki67','Cdkn1a','Krt8','Krt18','Krt7','Krt19','Sftpc','Lyz2','Sftpd','Sftpb','Abca3','Aqp5','Hopx','Vegfa'),]
newnames <- lapply(rownames(scmar),function(x) bquote(italic(.(x)))) 
colors <- rev(colorRampPalette(RColorBrewer::brewer.pal(10, "RdBu"))(256))
pheatmap(isr_counts[c('Mki67','Cdkn1a','Krt8','Krt18','Krt7','Krt19','Sftpc','Lyz2','Sftpd','Sftpb','Abca3','Aqp5','Hopx','Vegfa'),],
         annotation_col = dplyr::select(group_annot, c(Group,state)),
         annotation_colors = ann_colors,
         cluster_cols = F,
         cluster_rows = F,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         show_rownames = T,
         show_colnames = F,
         fontsize = 15,
         family = 'Arial',
         border_color = NA,
         #breaks = seq(from=-2,to=2,length.out=100),
         labels_row = as.expression(newnames),
         cellwidth = 20, cellheight = 20,
         gaps_col = c(7,14,21,26),
         scale = 'row',
         color = colors)
         
Atf=isr_counts[c('Atf3','Atf4','Atf5','Atf6'),]
newnames <- lapply(rownames(Atf),function(x) bquote(italic(.(x))))   
colors <- rev(colorRampPalette(RColorBrewer::brewer.pal(10, "RdBu"))(256))
pheatmap(isr_counts[c('Atf3','Atf4','Atf5','Atf6'),],
         annotation_col = dplyr::select(group_annot, c(Group,state)),
         annotation_colors = ann_colors,
         cluster_cols = F,
         cluster_rows = T,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         show_rownames = T,
         show_colnames = F,
         fontsize = 15,
         family = 'Arial',
         border_color = NA,
         #breaks = seq(from=-2,to=2,length.out=100),
         labels_row = as.expression(newnames),
         cellwidth = 20, cellheight = 20,
         gaps_col = c(7,14,21,26),
         scale = 'row',
         color = colors)
         
isr = read.table('/home/young/isr.txt')
isr = isr$V1 %>% unique()
isr_heatmap = isr_counts[isr,]
isr_heatmap = na.omit(isr_heatmap)
newnames <- lapply(isr,function(x) bquote(italic(.(x))))  
pheatmap(isr_heatmap,
         annotation_col = dplyr::select(group_annot, c(Group,state)),
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
         gaps_col = c(7,14,21,26),
         scale = 'row',
         color = colors)


