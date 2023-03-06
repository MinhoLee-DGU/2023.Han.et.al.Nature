##Remove batch effect
library(sva)
library(dplyr)
library(data.table)

raw_nd=read.table("./ndufs2_bulk/gene_read_counts_table_all_final.tsv", header=TRUE, stringsAsFactors=FALSE, row.names=1)
raw_isrib=read.table("./isrib_bulk/gene_read_counts_table_all_final.tsv", header=TRUE, stringsAsFactors=FALSE, row.names=1)
raw_nmn=read.table("./nmn_bulk/gene_read_counts_table_all_final.tsv", header=TRUE, stringsAsFactors=FALSE, row.names=1)

rawdata = cbind(raw_nd, raw_isrib, raw_nmn)
rawdata %>% head
rawdata = rawdata %>% as.matrix()
batch = c(rep(1,30),rep(2,64))

adjusted_counts <- ComBat_seq(rawdata, batch=batch, group=NULL, covar_mod=NULL)
write.table(adjusted_counts, "./ndufs2_isrib_nmn/adjusted_counts.tsv",sep = '\t', quote = F)


##Heatmap
library(pheatmap)
library(viridis)

counts = fread("./ndufs2_isrib_nmn/adjusted_counts.tsv",sep = '\t')
sample_info = read.csv('./ndufs2_isrib_nmn/sample_info.csv')

all_count=mat
library('edgeR')
class = colData$condition
genes=rownames(all_count)


y <- DGEList(counts=all_count, genes=genes, group=class)

y_filtered_norm <- calcNormFactors(y[filterByExpr(y),keep.lib.sizes = FALSE])

dge_design <- model.matrix(~Sex + class)
y_dsp <- estimateDisp(y_filtered_norm, dge_design, robust = TRUE)
y_qlf_test <- glmQLFit(y_dsp, dge_design, robust = TRUE)
y_qlf_test <- glmQLFTest(y_qlf_test)

group_annot <- sample_info %>% arrange(Group)
group_annot$Group = factor(group_annot$Group, levels = c('WT','NDUFS2 control','NDUFS2 control + ISRIB','NDUFS2 control + NMN','NDUFS2 control/NDI1','NDUFS2 cKO','NDUFS2 cKO + ISRIB','NDUFS2 cKO + NMN','NDUFS2 cKO/NDI1'))
group_annot$Sex = factor(group_annot$Sex, levels = c('F','M'))
group_annot = group_annot[order(group_annot$Group, group_annot$Sex),]

rownames(group_annot) = group_annot$Sample.name
group_annot = group_annot[,-1]

ann_colors = list(Group = c('WT' = '#000000', 'NDUFS2 control' = '#44000A', 'NDUFS2 control + ISRIB' = '#660010', 'NDUFS2 control + NMN' ='#B84A5B', 'NDUFS2 control/NDI1' = '#CF7C89',  
                            'NDUFS2 cKO' = '#005174', 'NDUFS2 cKO + ISRIB' = '#007AAE', 'NDUFS2 cKO + NMN' = '#60C5F1', 'NDUFS2 cKO/NDI1' = '#94D8F6'), Sex = c('F' = '#818181', 'M' = '#C0BEBF'))


isr_sum <- decideTests(y_qlf_test, p.value = 1)

isr_counts <- as.logical(abs(isr_sum)) %>%
  y_filtered_norm[.,] %>%
  edgeR::cpm(.) %>%
  as.data.frame() %>% 
  dplyr::select(rownames(group_annot)) 
  
Atf=isr_counts[c('Atf3','Atf4','Atf5','Atf6'),]
newnames <- lapply(rownames(Atf),function(x) bquote(italic(.(x))))   
colors <- rev(colorRampPalette(RColorBrewer::brewer.pal(10, "RdBu"))(256))
pheatmap(isr_counts[c('Atf3','Atf4','Atf5','Atf6'),],
         annotation_col = dplyr::select(group_annot, c(Group,Sex)),
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
         cellwidth = 10, cellheight = 20,
         gaps_col = c(6,27,35,46,53,66,75,86,94),
         scale = 'row',
         family = 'Arial',
         color = colors)
         
ndu=isr_counts[c('Ndufs2'),]
newnames <- lapply(rownames(ndu),function(x) bquote(italic(.(x))))   
colors <- rev(colorRampPalette(RColorBrewer::brewer.pal(10, "RdBu"))(256))
pheatmap(isr_counts[c('Ndufs2'),],
         annotation_col = dplyr::select(group_annot, c(Group,Sex)),
         annotation_colors = ann_colors,
         cluster_cols = F,
         cluster_rows = F,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         show_rownames = T,
         show_colnames = F,
         fontsize = 13,
         family = 'Arial',
         border_color = NA,
         #breaks = seq(from=-2,to=2,length.out=100),
         labels_row = as.expression(newnames),
         cellwidth = 10, cellheight = 40,
         gaps_col = c(6,27,35,46,53,66,75,86,94),
         scale = 'row',
         family = 'Arial',
         color = colors)

isr = read.table('./Bulk/isr.txt')
isr = isr$V1 %>% unique()
isr_heatmap = isr_count[isr,]
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
         gaps_col = c(6,27,35,46,53,66,75,86,94),
         scale = 'row',
         family = 'Arial',
         color = colors)
         
 


