# 2023.Han.et.al.Nature

R Notebooks for Han et al. 2023

All raw data (.fastq) from the sequencing analysis in this work are available at the NCBI BioProject with the following Accsession IDs:(scRNA-seq : PRJNA865889, Bulk-seq : PRJNA940730, PRJNA940746, PRJNA940973, PRJNA940986 and PRJNA940992)

scRNA : 
1. processing : Make GRCm39 reference genome with tdTomato inserted and remove doublet
2. whole_clusters : Annotate celltype for whole cells
3. subset_epithelial : Subclustering analysis with Epithelium
4. subset_Mesenchymal : Subclustering analysis with Mesencymal cells
5. integration_Hurskainen_Han.R : Epithelial cells from Hurskainen et al. data integrated with our scRNA data of epithelium
6. integration_Hurskainen_Han.py : UMAP embedding and RNA velocity prediction
7. integration_Negretti_Han.R : Epitehlial cells from Negretti et al. data integrated with our scRNA data of epithelium
8. integration_Negertti_Han.R :  UMAP embedding and RNA velocity prediction
9. other models and NDUFS2 cKO Epithelium : ISR enrichment scores across the epithelial cells within Strunz et al., Choi et al. and Kobayashi et al.

Bulk :
1. Ndufs2 cKO mice.R : Heatmap within NDUFS2 control, NDUFS2 control/NDI1, NDUFS2 cKO and NDUFS2 cKO/NDI1
2. Ndufs2 cKO, ISRIB and NMN.R : Heatmap within NDUFS2 cKO datasets, ISRIB datasets and NMN datasets
3. Ndufs2 cKO P35_P6 : Heatmap within P35 NDUFS2 cKO mice and P6 NDUFS2 cKO mice
4. Sdhd cKO and Ndufs2 cKO.R : Heatmap within SDHD control, SDHD cKO, NDUFS2 control and NDUFS2 cKO
5. P6_PCA_plot.R : PCA plot to visualize the clustering patterns of the samples whith in P6 NDUFS2 cKO datasets
6. P6_WT_cultured.R : Heatmap within incubated with piericidin A Cells and without piericidin A from P6 WT datasets
