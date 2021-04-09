#! /data/salomonis2/LabFiles/Frank-Li/citeseq/scanpy_new_env/bin/python3.6

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import scrublet as scr


# # scanpy workflow
# adata_all = sc.read_10x_mtx('../aat1699_DataS1/data/')
# adata_all.var['mt'] = adata_all.var_names.str.startswith('MT-')
# sc.pp.calculate_qc_metrics(adata_all,qc_vars=['mt'],percent_top=None,inplace=True,log1p=False)

# # first, only consider RCC patients and samples
# meta = pd.read_csv('../meta_sample.txt',sep='\t')
# need = meta.loc[(meta['Experiment'].str.startswith('RCC'))&(meta['Organ']=='Kidney'),:]
# need_sangerID = need['SangerID'].tolist()
# sangerID_to_patient = pd.Series(data=need['Experiment'].values,index=need['SangerID'].values).to_dict()
# sangerID_to_channel = pd.Series(data=[str(int(item)) for item in need['Channel10X'].values],index=need['SangerID'].values).to_dict()
# sangerID_to_status = pd.Series(data=need['TissueDiseaseState'].values,index=need['SangerID'].values).to_dict()


# adata_all.obs['SangerID'] = [item.split('___')[0] for item in adata_all.obs_names]
# adata = adata_all[adata_all.obs['SangerID'].isin(need_sangerID),:]
# adata.obs['patient'] = adata.obs['SangerID'].map(sangerID_to_patient).values
# adata.obs['channel'] = adata.obs['SangerID'].map(sangerID_to_channel).values
# adata.obs['status'] = adata.obs['SangerID'].map(sangerID_to_status).values


# # focus on certain patients
# adata = adata[(adata.obs['patient']=='RCC2')|(adata.obs['patient']=='RCC1'),:]



# # start normal analysis
# adata.obs.to_csv('./inspection_obs.txt',sep='\t')
# adata.var.to_csv('./inspection_var.txt',sep='\t')



# # run scrublet
# counts_matrix = adata.X.copy().toarray()
# print(counts_matrix.shape[0],counts_matrix.shape[1])
# scrub = scr.Scrublet(counts_matrix)
# doublet_scores,predicted_doublets = scrub.scrub_doublets(min_counts=1,min_cells=1)
# scrub.plot_histogram()
# plt.savefig('./scrub_check.pdf',bbox_inches='tight')
# plt.close()
# adata.obs['doublet_scores'] = doublet_scores
# adata.obs['predicted_doublets'] = np.logical_not(predicted_doublets)
# adata = adata[adata.obs['predicted_doublets'],:] 
# print(adata)

# sc.pp.filter_cells(adata,min_genes=200)
# sc.pp.filter_cells(adata,max_genes=8000)
# sc.pp.filter_genes(adata,min_cells=1)
# adata = adata[adata.obs['pct_counts_mt']<20,:]
# adata = adata[adata.obs['total_counts']<60000,:]

# sc.pp.normalize_total(adata,target_sum=1e4)
# sc.pp.log1p(adata)

# adata.raw = adata
# sc.pp.highly_variable_genes(adata,flavor='cell_ranger',n_top_genes=3000)
# adata = adata[:,adata.var['highly_variable']]
# sc.pp.regress_out(adata,['total_counts'])
# sc.pp.scale(adata,max_value=10)
# sc.tl.pca(adata)
# #sc.pp.neighbors(adata)
# sc.external.pp.bbknn(adata,batch_key='patient')
# sc.tl.leiden(adata)
# sc.tl.umap(adata)

# adata.write('./RCC12_after_umap_bbknn.h5ad')
# adata.raw.to_adata().write('./RCC12_to_cellxgene_bbknn.h5ad')


# plot interesting genes on UMAP
# adata = sc.read('./RCC2_after_umap.h5ad')
# sc.pl.umap(adata,color=['NDUFA5','NDUFAB1','NDUFAF3','NDUFS3','UQCRC1'],ncols=3)
# plt.savefig('./figures/ETC_genes.pdf',bbox_inches='tight')
# plt.close()
# sc.pl.umap(adata,color=['MRPL12','MRPL14','MRPL38','MRPL47','MRPS10','MRPS18A','MRPS24','MRPS36','MRPS5'],ncols=3)
# plt.savefig('./figures/MRP_genes.pdf',bbox_inches='tight')
# plt.close()
# sc.pl.umap(adata,color=['HLA-DMB','HLA-DOA','HLA-DPB1','HLA-DQA2','HLA-DQB2','HLA-DRA'],ncols=3)
# plt.savefig('./figures/MHCII_genes.pdf',bbox_inches='tight')
# plt.close()

# after annotation
adata = sc.read('./RCC12_after_umap_bbknn.h5ad')
annotation_map ={
    '0': 'non_cancer',
    '1': 'non_cancer',
    '2': 'non_cancer',
    '3': 'ccRCC',
    '4': 'non_cancer',
    '5': 'non_cancer',
    '6': 'non_cancer',
    '7': 'non_cancer',
    '8': 'non_cancer',
    '9': 'non_cancer',
    '10': 'non_cancer',
    '11': 'non_cancer',
    '12': 'non_cancer',
    '13': 'non_cancer',
    '14': 'ccRCC',
    '15': 'non_cancer',
    '16': 'non_cancer',
    '17': 'non_cancer',
    '18': 'non_cancer',
    '19': 'non_cancer',
    '20': 'non_cancer',
    '21': 'non_cancer',
    '22': 'non_cancer',
    '23': 'non_cancer',
    '24': 'non_cancer',
    '25': 'non_cancer',
    '26': 'non_cancer',
    '27': 'non_cancer',
}
adata.obs['anno'] = adata.obs['leiden'].astype('str').map(annotation_map)

# only focus on ccRCC
adata_c = adata[adata.obs['anno']=='ccRCC',:]
sc.tl.pca(adata_c)
#sc.pp.neighbors(adata_c)
sc.external.pp.bbknn(adata_c,batch_key='patient')
sc.tl.leiden(adata_c,resolution=0.5)
sc.tl.umap(adata_c)
#adata_c.raw.to_adata().write('./ccRCC_bbknn.h5ad')

# before DE, first get rid of uninteresting genes
cond1 = np.logical_not(adata_c.raw.var_names.str.startswith('MT-'))
cond2 = np.logical_not(adata_c.raw.var_names.isin(['MALAT1']))
is_interested = np.logical_and(cond1,cond2)
adata_c.raw = adata_c.raw[:,is_interested].to_adata()

# different gene analysis, automatic genes
sc.tl.rank_genes_groups(adata_c,groupby='leiden')
sc.pl.rank_genes_groups(adata_c,save='.pdf')
sc.pl.rank_genes_groups_heatmap(adata_c,n_genes=8,swap_axes=True,save='.pdf')

# different gene analysis, interesting genes
genes = ['NDUFA5','NDUFAB1','NDUFAF3','NDUFS3','UQCRC1','MRPL12','MRPL14','MRPL38','MRPL47',
        'MRPS10','MRPS18A','MRPS24','MRPS36','MRPS5','HLA-DMB','HLA-DOA','HLA-DPB1','HLA-DQA2','HLA-DQB2','HLA-DRA']
sc.pl.heatmap(adata_c,var_names=genes,groupby='leiden',swap_axes=True)
plt.savefig('./figures/21gene_de.pdf',bbox_inches='tight')
plt.close()

# # different state analysis, automatic
# for cluster in adata.obs['anno'].astype('category').cat.categories:
#     adata_c = adata[adata.obs['anno']==cluster,:]
#     try:
#         sc.tl.rank_genes_groups(adata_c,groupby='status')
#     except:
#         print("cluster {} is in one status\n".format(cluster))
#         continue
#     else:
#         sc.pl.rank_genes_groups(adata_c)
#         plt.savefig('./figures/{0}_elbow.pdf'.format(cluster),bbox_inches='tight')
#         plt.close()
#         sc.pl.rank_genes_groups_heatmap(adata_c,n_genes=15,swap_axes=True)
#         plt.savefig('./figures/{0}_heatmap.pdf'.format(cluster),bbox_inches='tight')
#         plt.close()

# # differential state analysis, interesting genes
# for cluster in adata.obs['anno'].astype('category').cat.categories:
#     adata_c = adata[adata.obs['anno']==cluster,:]
#     try:
#         sc.pl.heatmap(adata_c,var_names=genes,groupby='status',swap_axes=True)
#     except:
#         print("cluster {} is in one status\n".format(cluster))
#         continue
#     else:
#         plt.savefig('./figures/{}_21gene_heatmap.pdf'.format(cluster),bbox_inches='tight')
#         plt.close()

# # only need ccRCC and 23 genes, and clustering on that
# adata_i = adata.raw.to_adata()[adata.obs['anno']=='ccRCC',genes]
# df = pd.DataFrame(data=adata_i.X.copy().toarray().T,index=adata_i.var_names)
# sns.clustermap(data=df,metric='euclidean',cmap='viridis',row_cluster=False)
# plt.savefig('./figures/cluserheatmap.pdf',bbox_inches='tight')
# plt.close()


#adata.raw.to_adata().write('./RCC2_to_cellxgene.h5ad')

