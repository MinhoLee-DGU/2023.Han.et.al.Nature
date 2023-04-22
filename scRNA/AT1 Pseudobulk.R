library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(DESeq2)
library(tidyverse)
library(ggrepel)
library(edgeR)
library(ggplot2)
library(ggrepel)

epi = readRDS('../epi.rds')


epi@meta.data$samples = paste0(epi@meta.data$state, '_', epi@meta.data$orig.ident)
cts <- AggregateExpression(epi, 
                    group.by = c("chart", "samples"),
                    assays = 'RNA',
                    slot = "counts",
                    return.seurat = FALSE)

cts <- cts$RNA
cts.t <- t(cts)
cts.t <- as.data.frame(cts.t)
splitRows <- gsub('_.*', '', rownames(cts.t))

cts.split <- split.data.frame(cts.t,
                 f = factor(splitRows))

cts.split.modified <- lapply(cts.split, function(x){
  rownames(x) <- gsub('.*_(.*)', '\\1', rownames(x))
  t(x) 
})

counts_at1 <- cts.split.modified$'AT1'


#generate colData
colData <- data.frame(samples = colnames(counts_at1))

colData$condition = rep(c('NDUFS2 cKO','NDUFS2 control'),c(4,4))
colData$condition = factor(colData$condition, levels = c('NDUFS2 control','NDUFS2 cKO'))

# perform DESeq2 
# Create DESeq2 object   
dds <- DESeqDataSetFromMatrix(countData =counts_at1,
                       colData = colData,
                       design = ~ condition)

# filter
keep <- rowSums(counts(dds)) >=10
dds <- dds[keep,]

# run DESeq2
dds <- DESeq(dds)
rld = dds %>% rlog()
DESeq2::plotPCA(rld,intgroup=c("condition"))

colData(rld)$label = rld@colData %>% rownames

DESeq2::plotPCA(rld, colnames(colData(rld)), returnData = T) %>% 
  ggplot(aes(PC1, PC2, color = condition, label = label)) +
  geom_label_repel(aes(PC1, PC2, label = colData(rld)$label)) +
  geom_point(size = 4) + 
  xlab(paste0("PC1: 50% variance")) +
  ylab(paste0("PC2: 21% variance")) + 
  theme_bw()


#using edgeR
class = colData$condition
genes=rownames(counts_at1)

y <- DGEList(counts=counts_at1, genes=genes, group=class)

y_filtered_norm <- calcNormFactors(y[filterByExpr(y),keep.lib.sizes = FALSE])

dge_design <- model.matrix(~class)
y_dsp <- estimateDisp(y_filtered_norm, dge_design, robust = TRUE)
y_qlf_test <- glmQLFit(y_dsp, dge_design, robust = TRUE)
y_qlf_test <- glmQLFTest(y_qlf_test)

y_top_tags <- topTags(y_qlf_test, n = "inf")
y_top_tags
DE_genes = y_top_tags %>% as.data.frame()

#volcanoplot
DE_genes$group = 'No.sig'
DE_genes$group[DE_genes$FDR < 0.05 & DE_genes$logFC > 1.5] = 'NDUFS2 cKO'
DE_genes$group[DE_genes$FDR < 0.05 & DE_genes$logFC < -1.5] = 'NDUFS2 control'

cols = c("#00BFC4",'#BFBFBF','#F8766D')
names(cols) = c("NDUFS2 cKO","No.sig","NDUFS2 control")

DE_genes$label = 'No'
DE_genes$label[DE_genes$genes == 'Igfbp2'] = 'label'
DE_genes$label[DE_genes$genes == 'Ndufs2'] = 'label'
DE_genes$label[DE_genes$genes == 'Eif4ebp1'] = 'label'
DE_genes$label[DE_genes$genes == 'Atf5'] = 'label'
DE_genes$label[DE_genes$genes == 'Phgdh'] = 'label'
DE_genes$label[DE_genes$genes == 'Krt18'] = 'label'
DE_genes$label[DE_genes$genes == 'Krt8'] = 'label'
DE_genes$label[DE_genes$genes == 'Cdkn1a'] = 'label'

ggplot(DE_genes, aes(x = logFC, y = -log10(FDR), color = group)) +
scale_colour_manual(values = cols) +
ggtitle(label = "Volcano Plot (AT1)", subtitle = "NDUFS2 control vs NDUFS2 cKO") +
geom_point(size = 2.5, alpha = 0.7, na.rm = T) + geom_text_repel(aes(x = logFC, y = -log10(FDR), label = ifelse(label == 'label', genes,"")),alpha=2,family ='Arial', box.padding = 0.5, max.overlaps = Inf, size = 4) +
theme_bw(base_size = 14) + 
theme(legend.position = "right") + 
xlab("logFC") + 
ylab(expression(-log[10]("FDR"))) +
geom_vline(xintercept = 1.5, colour="#8A8A8A", linetype="dashed") + 
geom_vline(xintercept = -1.5, colour="#8A8A8A", linetype="dashed")+ 
scale_y_continuous(trans = "log1p") + ylim(c(0,6)) + xlim(c(-6,6)) +
theme(plot.title = element_text(size = 20,family = 'Arial'), legend.text = element_text(size = 15,family = 'Arial')) + 
theme(axis.title.x = element_text(size = 15,family = 'Arial'), axis.title.y = element_text(size = 15,family = 'Arial')) + 
theme(axis.text.x = element_text(size = 15, family = 'Arial'), axis.text.y = element_text(size = 15,family = 'Arial'))



