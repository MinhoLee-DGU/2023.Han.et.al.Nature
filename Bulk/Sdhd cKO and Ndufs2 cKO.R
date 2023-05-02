##Remove batch effect
library(sva)
library(dplyr)
library(data.table)
library(edgeR)
library(pheatmap)
library(viridis)

raw_nd=read.table("./ndufs2_bulk/gene_read_counts_table_all_final.tsv", header=TRUE, stringsAsFactors=FALSE, row.names=1)
raw_sd=read.table("./sdhd_bulk/gene_read_counts_table_all_final.tsv", header=TRUE, stringsAsFactors=FALSE, row.names=1)

rawdata = cbind(raw_nd, raw_sd)
rawdata %>% head
rawdata = rawdata %>% as.matrix()
batch = c(rep(1,30),rep(2,13))

adjusted_counts = ComBat_seq(rawdata, batch=batch, group=NULL, covar_mod=NULL)
write.table(adjusted_counts, "./nd_sd/adjusted_counts.tsv",sep = '\t', quote = F)



##Heatmap
counts = fread("./nd_sd/adjusted_counts.tsv",sep = '\t')
mat = counts %>% as.data.frame()
rownames(mat) = mat$V1
mat = mat[,-1]

#only 'Sdhd Control','Sdhd cKO','Ndufs2 Control' and 'Ndufs2 cKO' sample
sample_info = read.csv('./nd_sd/sample_info.csv')

mat = mat[,sample_info$Sample.name]

class = factor(sample_info$Group, levels = c('Sdhd Control','Sdhd cKO','Ndufs2 Control','Ndufs2 cKO'))
genes=rownames(mat)
Sex = factor(sample_info$sex, levels=c('F','M'))

y = DGEList(counts=mat, genes=genes, group=class)
y_filtered_norm = calcNormFactors(y[filterByExpr(y),keep.lib.sizes = FALSE])

dge_design = model.matrix(~Sex + class)
y_dsp = estimateDisp(y_filtered_norm, dge_design, robust = TRUE)
y_qlf_test = glmQLFit(y_dsp, dge_design, robust = TRUE)
y_qlf_test = glmQLFTest(y_qlf_test)

group_annot = sample_info %>% arrange(Group)
group_annot$Group = factor(group_annot$Group, levels = c('Sdhd Control','Sdhd cKO','Ndufs2 Control','Ndufs2 cKO'))
group_annot$sex = factor(group_annot$sex, levels = c('F','M'))
group_annot = group_annot[order(group_annot$Group, group_annot$sex),]

ann_colors = list(Group = c('Sdhd Control' = '#7CAE00', 'Sdhd cKO' = '#C77CFF', 'Ndufs2 Control' = '#00BFC4',  'Ndufs2 cKO' = '#F8766D'), sex = c('F' = '#818181', 'M' = '#C0BEBF'))

isr_sum <- decideTests(y_qlf_test, p.value = 1)
rownames(group_annot) = group_annot$Sample.name
group_annot = group_annot[,-1]

isr_counts <- as.logical(abs(isr_sum)) %>%
  y_filtered_norm[.,] %>%
  edgeR::cpm(.) %>%
  as.data.frame() %>% 
  dplyr::select(rownames(group_annot))

#### Convert rownames of isr_counts to Gene Symbol !!!! ####
isr_heatmap = isr_counts # isr_counts rownames changed to the genetic symbols 

Atf=isr_heatmap[c('Atf3','Atf4','Atf5','Atf6'),]
newnames <- lapply(rownames(Atf),function(x) bquote(italic(.(x))))   

colors <- rev(colorRampPalette(RColorBrewer::brewer.pal(10, "RdBu"))(256))
pheatmap(isr_heatmap[c('Atf3','Atf4','Atf5','Atf6'),],
         annotation_col = dplyr::select(group_annot, c(Group,sex)),
         annotation_colors = ann_colors,
         cluster_cols = F,
         cluster_rows = T,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         show_rownames = T,
         show_colnames = F,
         fontsize = 28,
         family = 'Arial',
         border_color = NA,
         #breaks = seq(from=-2,to=2,length.out=100),
         labels_row = as.expression(newnames),
         cellwidth = 20, cellheight = 40,
         gaps_col = c(6,13,21,28),
         scale = 'row',
         color = colors)


isr = read.table('./isr.txt')
isr = isr$V1 %>% unique()
isr_heatmap = isr_counts[isr,]
isr_heatmap = na.omit(isr_heatmap)

colors <- rev(colorRampPalette(RColorBrewer::brewer.pal(10, "RdBu"))(256))
newnames <- lapply(isr,function(x) bquote(italic(.(x))))  
pheatmap(isr_heatmap,
         annotation_col = dplyr::select(group_annot, c(Group,sex)),
         annotation_colors = ann_colors,
         cluster_cols = F,
         cluster_rows = T,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         show_rownames = F,
         show_colnames = F,
         fontsize = 14,
         border_color = NA,
         fontsize_row = 6,
         family = 'Arial',
         #breaks = seq(from=-2,to=2,length.out=100),
         #labels_row = as.expression(newnames),
         #cellwidth = 20, cellheight = 40,
         gaps_col = c(6,13,21,28),
         scale = 'row',
         family = 'Arial',
         color = colors)
 




