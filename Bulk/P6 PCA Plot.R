library(DESeq2)
library(tidyverse)
library(ggrepel)
library(edgeR)


##PCA plot
dir = './p6_bulk/'
sample_info = read.csv('./p6_bulk/sample_info_p6.csv')
sample_info$Group[sample_info$Group == 'Control'] = 'P6 NDUFS2 control'
sample_info$Group[sample_info$Group == 'KO'] = 'P6 NDUFS2 cKO'
pca = sample_info
rownames(pca) = pca$Sample.name
pca = pca[c('P6_1','P6_10','P6_11','P6_12','P6_2','P6_3','P6_4','P6_5','P6_6','P6_7','P6_8','P6_9'),]

sampleFiles = grep("A*_gene.tsv",list.files(dir),value=TRUE)
sampleCondition = sub("(*).gene.*","\\1",sampleFiles)
sampleTable = data.frame(sampleName = sampleFiles,  fileName = sampleFiles,  condition = sampleCondition)

sampleTable$condition = factor(sampleTable$condition) 
ddsHTSeq = DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,  directory = dir, design= ~ condition)
keep = rowSums(counts(ddsHTSeq)) >= 10
dds = ddsHTSeq[keep,]
dds$condition = factor(pca$Group, levels = c('P6 NDUFS2 control','P6 NDUFS2 cKO'))
dds = DESeq(dds)
rld = dds %>% rlog()

colData(rld)$sex = factor(sample_info$Sex, levels = c('F', 'M'))
colData(rld)$label = sample_info$Sample.name
Sex = factor(sample_info$Sex, levels = c('F','M'))


p=DESeq2::plotPCA(rld, colnames(colData(rld)), returnData = T) %>% 
  ggplot(aes(PC1, PC2,shape = Sex),colour="black")  +
  geom_point(aes(color = condition), size=10) +
  geom_point(shape = Sex,size = 10,colour = "black") +
  scale_color_manual(values = c('#660010','#007AAE'), name = '') +
  xlab(paste0("PC1: 38% variance")) +
  ylab(paste0("PC2: 26% variance")) + 
  theme_bw() +
  theme(plot.title = element_text(size=15,family = 'Arial'), legend.text = element_text(size=15,family = 'Arial')) + 
  theme(axis.title.x = element_text(size=15,family = 'Arial'), axis.title.y = element_text(size=15,family = 'Arial')) + 
  theme(axis.text.y = element_text(size=15,family = 'Arial')) + theme(axis.text.x = element_text(size=15,family = 'Arial'))

p + theme(legend.position='bottom') + guides(color = guide_legend(nrow = 2, byrow = T, override.aes = list(size=10)))







