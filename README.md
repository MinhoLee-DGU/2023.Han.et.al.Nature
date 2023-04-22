# 2023.Han.et.al.Nature

R Notebooks for Han et al. 2023

All raw data (.fastq) from the sequencing analysis in this work are available at the NCBI BioProject with the following Accsession IDs:(scRNA-seq : PRJNA865889, Bulk-seq : PRJNA940730, PRJNA940746, PRJNA940973, PRJNA940986 and PRJNA940992)

scRNA : 
1. Processing : Make GRCm39 reference genome with tdTomato inserted and remove doublet
2. Whole clusters : Annotate celltype for whole cells
3. Subset epithelial : Subclustering analysis with Epithelium
4. Subset Mesenchymal : Subclustering analysis with Mesencymal cells
5. Subset AT1 cells : subclustering analysis with AT1 cells
6. Integration Hurskainen_Han.R : Epithelial cells from Hurskainen et al. data integrated with our scRNA data of epithelium
7. Integration Hurskainen_Han.py : UMAP embedding and RNA velocity prediction
8. Integration Negretti_Han.R : Epitehlial cells from Negretti et al. data integrated with our scRNA data of epithelium
9. Integration Negertti_Han.R :  UMAP embedding and RNA velocity prediction
10. Other models and NDUFS2 cKO Epithelium : ISR enrichment scores across the epithelial cells within Strunz et al., Choi et al. and Kobayashi et al.

Bulk :
1. NDUFS2 cKO mice.R : Heatmap within NDUFS2 control, NDUFS2 control/NDI1, NDUFS2 cKO and NDUFS2 cKO/NDI1
2. NDUFS2 cKO, ISRIB and NMN.R : Heatmap within NDUFS2 cKO datasets, ISRIB datasets and NMN datasets
3. NDUFS2 cKO P35 P6 : Heatmap within P35 NDUFS2 cKO mice and P6 NDUFS2 cKO mice
4. SDHD cKO and NDUFS2 cKO.R : Heatmap within SDHD control, SDHD cKO, NDUFS2 control and NDUFS2 cKO
5. P6 PCA plot.R : PCA plot to visualize the clustering patterns of the samples whith in P6 NDUFS2 cKO datasets
6. P6 WT cultured.R : Heatmap within incubated with piericidin A Cells and without piericidin A from P6 WT datasets
