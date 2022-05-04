##### scanpy in python
import anndata
import matplotlib as plt
import scvelo as scv
import scanpy as sc
import cellrank as cr
import numpy as np
import pandas as pd
import anndata as ad
import bbknn

sample_obs = pd.read_csv("../integ_Neg_Han/cellID_obs.csv")
cell_clusters = pd.read_csv("../integ_Neg_Han/clusters.csv")
Age = pd.read_csv('../integ_Neg_Han/Age.csv')
project = pd.read_csv('../integ_Neg_Han/project.csv')
cell_colors = pd.read_csv('../integ_Neg_Han/colors.csv')


# load loom files from current work
sample1 = anndata.read_loom("/home/young/scRNA_study/cell_ranger/Mus_musculus/only_tdt/SHH1/velocyto/SHH1.loom")
sample2 = anndata.read_loom("/home/young/scRNA_study/cell_ranger/Mus_musculus/only_tdt/SHH2/velocyto/SHH2.loom")
sample3 = anndata.read_loom("/home/young/scRNA_study/cell_ranger/Mus_musculus/only_tdt/SHH3/velocyto/SHH3.loom")
sample4 = anndata.read_loom("/home/young/scRNA_study/cell_ranger/Mus_musculus/only_tdt/SHH4/velocyto/SHH4.loom")
sample5 = anndata.read_loom("/home/young/scRNA_study/cell_ranger/Mus_musculus/only_tdt/SHH5/velocyto/SHH5.loom")
sample6 = anndata.read_loom("/home/young/scRNA_study/cell_ranger/Mus_musculus/only_tdt/SHH6/velocyto/SHH6.loom")
sample7 = anndata.read_loom("/home/young/scRNA_study/cell_ranger/Mus_musculus/only_tdt/SHH7/velocyto/SHH7.loom")
sample8 = anndata.read_loom("/home/young/scRNA_study/cell_ranger/Mus_musculus/only_tdt/SHH8/velocyto/SHH8.loom")

# load loom files from Negretti et al.
sample9 = anndata.read_loom('/home/young/scRNA_study/21_ML/ref_ml/ref_biorxiv/RAW/P3/3890-JS-2/velocyto/3890-JS-2.loom')
sample10 = anndata.read_loom('/home/young/scRNA_study/21_ML/ref_ml/ref_biorxiv/RAW/P5/SRR13481279/3886-JB-1/velocyto/3886-JB-1.loom')
sample11 = anndata.read_loom('/home/young/scRNA_study/21_ML/ref_ml/ref_biorxiv/RAW/P5/SRR13481280/3886-JB-3/velocyto/3886-JB-3.loom')
sample12 = anndata.read_loom('/home/young/scRNA_study/21_ML/ref_ml/ref_biorxiv/RAW/P7/SRR12990750/3890-JS-3/velocyto/3890-JS-3.loom')
sample13 = anndata.read_loom('/home/young/scRNA_study/21_ML/ref_ml/ref_biorxiv/RAW/P7/SRR13481281/3890-JS-9/velocyto/3890-JS-9.loom')
sample14 = anndata.read_loom('/home/young/scRNA_study/21_ML/ref_ml/ref_biorxiv/RAW/P14/SRR12990752/3890-JS-4/velocyto/3890-JS-4.loom')
sample15 = anndata.read_loom('/home/young/scRNA_study/21_ML/ref_ml/ref_biorxiv/RAW/P14/SRR12990751/3890-JS-10/velocyto/3890-JS-10.loom')
sample16 = anndata.read_loom('/home/young/scRNA_study/21_ML/ref_ml/ref_biorxiv/RAW/P0/SRR12990747/3890-JS-1//velocyto/3890-JS-1.loom')
sample17 = anndata.read_loom('/home/young/scRNA_study/21_ML/ref_ml/ref_biorxiv/RAW/P0/SRR12990748/3890-JS-6/velocyto/3890-JS-6.loom')
sample18 = anndata.read_loom('/home/young/scRNA_study/21_ML/ref_ml/ref_biorxiv/RAW/P7/SRR12990749/3890-JS-13//velocyto/3890-JS-13.loom')

sample1.var_names_make_unique()
sample2.var_names_make_unique()
sample3.var_names_make_unique()
sample4.var_names_make_unique()
sample5.var_names_make_unique()
sample6.var_names_make_unique()
sample7.var_names_make_unique()
sample8.var_names_make_unique()
sample9.var_names_make_unique()
sample10.var_names_make_unique()
sample11.var_names_make_unique()
sample12.var_names_make_unique()
sample13.var_names_make_unique()
sample14.var_names_make_unique()
sample15.var_names_make_unique()
sample16.var_names_make_unique()
sample17.var_names_make_unique()
sample18.var_names_make_unique()

cellID_obs_sample1 = sample_obs[sample_obs["x"].str.contains("SHH1")]
cellID_obs_sample2 = sample_obs[sample_obs["x"].str.contains("SHH2")]
cellID_obs_sample3 = sample_obs[sample_obs["x"].str.contains("SHH3")]
cellID_obs_sample4 = sample_obs[sample_obs["x"].str.contains("SHH4")]
cellID_obs_sample5 = sample_obs[sample_obs["x"].str.contains("SHH5")]
cellID_obs_sample6 = sample_obs[sample_obs["x"].str.contains("SHH6")]
cellID_obs_sample7 = sample_obs[sample_obs["x"].str.contains("SHH7")]
cellID_obs_sample8 = sample_obs[sample_obs["x"].str.contains("SHH8")]
cellID_obs_sample9 = sample_obs[sample_obs["x"].str.contains("3890-JS-2")]
cellID_obs_sample10 = sample_obs[sample_obs["x"].str.contains("3886-JB-1")]
cellID_obs_sample11 = sample_obs[sample_obs["x"].str.contains("3886-JB-3")]
cellID_obs_sample12 = sample_obs[sample_obs["x"].str.contains("3890-JS-3")]
cellID_obs_sample13 = sample_obs[sample_obs["x"].str.contains("3890-JS-9")]
cellID_obs_sample14 = sample_obs[sample_obs["x"].str.contains("3890-JS-4")]
cellID_obs_sample15 = sample_obs[sample_obs["x"].str.contains("3890-JS-10")]
cellID_obs_sample16 = sample_obs[sample_obs["x"].str.contains("3890-JS-1")]
cellID_obs_sample17 = sample_obs[sample_obs["x"].str.contains("3890-JS-6")]
cellID_obs_sample18 = sample_obs[sample_obs["x"].str.contains("3890-JS-13")]

sample1 = sample1[np.isin(sample1.obs.index, cellID_obs_sample1)]
sample2 = sample2[np.isin(sample2.obs.index, cellID_obs_sample2)]
sample3 = sample3[np.isin(sample3.obs.index, cellID_obs_sample3)]
sample4 = sample4[np.isin(sample4.obs.index, cellID_obs_sample4)]
sample5 = sample5[np.isin(sample5.obs.index, cellID_obs_sample5)]
sample6 = sample6[np.isin(sample6.obs.index, cellID_obs_sample6)]
sample7 = sample7[np.isin(sample7.obs.index, cellID_obs_sample7)]
sample8 = sample8[np.isin(sample8.obs.index, cellID_obs_sample8)]
sample9 = sample9[np.isin(sample9.obs.index, cellID_obs_sample9)]
sample10 = sample10[np.isin(sample10.obs.index, cellID_obs_sample10)]
sample11 = sample11[np.isin(sample11.obs.index, cellID_obs_sample11)]
sample12 = sample12[np.isin(sample12.obs.index, cellID_obs_sample12)]
sample13 = sample13[np.isin(sample13.obs.index, cellID_obs_sample13)]
sample14 = sample14[np.isin(sample14.obs.index, cellID_obs_sample14)]
sample15 = sample15[np.isin(sample15.obs.index, cellID_obs_sample15)]
sample16 = sample16[np.isin(sample16.obs.index, cellID_obs_sample16)]
sample17 = sample17[np.isin(sample17.obs.index, cellID_obs_sample17)]
sample18 = sample18[np.isin(sample18.obs.index, cellID_obs_sample18)]

# merge sample
sample = sample1.concatenate(sample2, sample3, sample4, sample5, sample6, sample7, sample8, sample9, sample10, sample11, sample12, sample13, sample14, sample15, sample16, sample17, sample18)

sample_index = pd.DataFrame(sample.obs.index)
sample = sample_index.rename(columns = {0:'CellID'})

for i in range(0, 12394):
	sample.loc[i,'CellID'] = sample.loc[i,'CellID'][:-2]
  
for i in range(12394, 21129):
	sample.loc[i,'CellID'] = sample.loc[i,'CellID'][:-3]

# input cells info.
sample = sample.rename(columns = {'CellID':'Cell ID'})
 
cell_clusters = cell_clusters.rename(columns = {'Unnamed: 0':'Cell ID'})
clusters_ordered = sample.merge(cell_clusters, on = "Cell ID")
clusters_ordered = clusters_ordered.iloc[:,1:]

Age = Age.rename(columns = {'Unnamed: 0':'Cell ID'})
Age_ordered = sample.merge(Age, on = "Cell ID")
Age_ordered = Age_ordered.iloc[:,1:]

project = project.rename(columns = {'Unnamed: 0':'Cell ID'})
project_ordered = sample.merge(project, on = "Cell ID")
project_ordered = project_ordered.iloc[:,1:]

cell_colors = cell_colors.rename(columns = {'Unnamed: 0':'Cell ID'})
colors_ordered = sample.merge(cell_colors, on = "Cell ID")
colors_ordered = colors_ordered.iloc[:,1:]


sample = sample1.concatenate(sample2, sample3, sample4, sample5, sample6, sample7, sample8, sample9, sample10, sample11, sample12, sample13, sample14, sample15, sample16, sample17, sample18)
sample.uns['Cluster_colors'] = colors_ordered.values
sample.obs['Clusters'] = clusters_ordered.values
sample.obs['Age'] = Age_ordered.values
sample.obs['project'] = project_ordered.values


scv.pp.filter_and_normalize(sample)
scv.pp.moments(sample)
scv.tl.velocity(sample, mode = "stochastic")
scv.tl.velocity_graph(sample)
sc.tl.pca(sample)

# blancing batch effect
sc.external.pp.bbknn(sample, batch_key='Age')
sc.tl.umap(sample)

# visualize RNA velo city
color = sample.uns['Cluster_colors']
scv.pl.velocity_embedding_stream(sample, basis="umap", color=color, legend_loc = 'on data', legend_fontsize=12)




