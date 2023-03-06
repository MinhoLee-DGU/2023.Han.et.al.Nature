##Remove batch effect
library(sva)
library(dplyr)
library(data.table)
library(DESeq2)
library(tidyverse)
library(ggrepel)
library(pheatmap)
library(viridis)

raw_nd=read.table("./ndufs2_bulk/gene_read_counts_table_all_final.tsv", header=TRUE, stringsAsFactors=FALSE, row.names=1)
raw_p6=read.table("./p6_bulk/gene_read_counts_table_all_final.tsv", header=TRUE, stringsAsFactors=FALSE, row.names=1)

rawdata = cbind(raw_nd, raw_p6)
rawdata %>% head
rawdata = rawdata %>% as.matrix()
batch = c(rep(1,30),rep(2,12))

adjusted_counts <- ComBat_seq(rawdata, batch=batch, group=NULL, covar_mod=NULL)
write.table(adjusted_counts, "./nd_p6/adjusted_counts.tsv",sep = '\t', quote = F)


##PCA plot
counts = fread("./nd_p6/adjusted_counts.tsv" ,sep = '\t')

mat = counts
mat = mat %>% as.data.frame()
rownames(mat) = counts$GeneID
colData <- as.data.frame(colnames(mat))

sample_info = read.csv('./nd_p6/sample_info_nd_p6.csv')
sample_info = filter(sample_info, Group == 'ko_nd' | Group == 'Control' | Group == 'KO' | Group == 'con_nd')

sample_info$Group[sample_info$Group == 'ko_nd'] = 'P35 NDUFS2 cKO'
sample_info$Group[sample_info$Group == 'Control'] = 'P6 NDUFS2 control'
sample_info$Group[sample_info$Group == 'KO'] = 'P6 NDUFS2 cKO'
sample_info$Group[sample_info$Group == 'con_nd'] = 'P35 NDUFS2 control'

colData <- as.data.frame(colnames(mat[,sample_info$Sample.name]))

colData$condition = factor(sample_info$Group, levels = c('P6 NDUFS2 control','P35 NDUFS2 control','P6 NDUFS2 cKO','P35 NDUFS2 cKO'))
colnames(colData) <- c("colData", "condition")
ddsMat = DESeqDataSetFromMatrix(countData = mat[,sample_info$Sample.name], colData = colData, design = ~condition)

keep <- rowSums(counts(ddsMat)) >= 10
dds <- ddsMat[keep,]
dds <- DESeq(dds)
rld = dds %>% rlog()
DESeq2::plotPCA(rld,intgroup=c("condition"))


Sex = factor(sample_info$Sex, levels = c('F','M'))
colData(rld)$label = rld@colData$condition


  p=DESeq2::plotPCA(rld, colnames(colData(rld)), returnData = T) %>% 
  ggplot(aes(PC1, PC2,shape = Sex),colour="black")  +
  geom_point(aes(color = condition), size=10) +
  geom_point(shape = Sex,size = 10,colour = "black") +
  scale_color_manual(values = c('#660010','#CF7C89','#007AAE','#94D8F6'), name = '') +
  xlab(paste0("PC1: 43% variance")) +
  ylab(paste0("PC2: 19% variance")) + 
  theme_bw() +
  theme(plot.title = element_text(size=15,family = 'Arial'), legend.text = element_text(size=15,family = 'Arial')) + 
  theme(axis.title.x = element_text(size=15,family = 'Arial'), axis.title.y = element_text(size=15,family = 'Arial')) + 
  theme(axis.text.y = element_text(size=15,family = 'Arial')) + theme(axis.text.x = element_text(size=15,family = 'Arial'))

p + theme(legend.position='bottom') + guides(color = guide_legend(nrow = 2, byrow = T, override.aes = list(size=10)))

##Heatmap
all_count=mat[,sample_info$Sample.name]

library('edgeR')
class = colData$condition
genes=mat$V1

rownames(all_count) = genes


y <- DGEList(counts=all_count, genes=genes, group=class)

y_filtered_norm <- calcNormFactors(y[filterByExpr(y),keep.lib.sizes = FALSE])

dge_design <- model.matrix(~Sex + class)
y_dsp <- estimateDisp(y_filtered_norm, dge_design, robust = TRUE)
y_qlf_test <- glmQLFit(y_dsp, dge_design, robust = TRUE)
y_qlf_test <- glmQLFTest(y_qlf_test)

sample_info$Group = factor(sample_info$Group, levels = c('P6 NDUFS2 control','P35 NDUFS2 control','P6 NDUFS2 cKO','P35 NDUFS2 cKO'))
sample_info$Sex = factor(sample_info$Sex, levels = c('F','M'))

group_annot = sample_info
group_annot = group_annot[order(group_annot$Group, group_annot$Sex),]

rownames(group_annot) = group_annot$Sample.name
group_annot = group_annot[,-1]

isr_sum <- decideTests(y_qlf_test, p.value = 1)

isr_counts <- as.logical(abs(isr_sum)) %>%
  y_filtered_norm[.,] %>%
  edgeR::cpm(.) %>%
  as.data.frame() %>% 
  dplyr::select(rownames(group_annot))
 
ann_colors = list(Group = c('P6 NDUFS2 control' = '#660010', 'P35 NDUFS2 control' = '#CF7C89', 'P6 NDUFS2 cKO' = '#007AAE', 'P35 NDUFS2 cKO' = '#94D8F6'),Sex = c('F' = '#818181', 'M' = '#C0BEBF'))

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
         cellwidth = 20, cellheight = 30,
         gaps_col = c(6,14,20,27),
         scale = 'row',
         family = 'Arial',
         color = colors)

isr = read.table('./Bulk/isr.txt')
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
         fontsize = 15,
         border_color = NA,
         fontsize_row = 15,
         family = 'Arial',
         #breaks = seq(from=-2,to=2,length.out=100),
         #labels_row = as.expression(newnames),
         #cellwidth = 20, cellheight = 40,
         gaps_col = c(6,14,20,27),
         scale = 'row',
         family = 'Arial',
         color = colors)

