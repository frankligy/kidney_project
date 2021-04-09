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


# # focus on single patient
# adata = adata[adata.obs['patient']=='RCC2',:]

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
# sc.pp.filter_genes(adata,min_cells=1)
# adata = adata[adata.obs['pct_counts_mt']<20,:]
# adata = adata[adata.obs['total_counts']<40000,:]

# sc.pp.normalize_total(adata,target_sum=1e4)
# sc.pp.log1p(adata)

# adata.raw = adata
# sc.pp.highly_variable_genes(adata,flavor='cell_ranger',n_top_genes=3000)
# adata = adata[:,adata.var['highly_variable']]
# sc.pp.regress_out(adata,['total_counts'])
# sc.pp.scale(adata,max_value=10)
# sc.tl.pca(adata)
# sc.pp.neighbors(adata)
# sc.tl.leiden(adata)
# sc.tl.umap(adata)

# adata.write('./RCC2_after_umap.h5ad')
# adata.raw.to_adata().write('./RCC2_to_cellxgene.h5ad')


# plot interesting genes on UMAP
adata = sc.read('./RCC2_after_umap.h5ad')
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
annotation_map ={
    '0': 'CD4_T',
    '1': 'CD8_T',
    '2': 'Macrophage',
    '3': 'ccRCC',
    '4': 'proximal',
    '5': 'proximal',
    '6': 'CD4_T',
    '7': 'Macrophage',
    '8': 'glomerular_endothelial',
    '9': 'Tregs',
    '10': 'NK',
    '11': 'glomerular_endothelial',
    '12': 'ccRCC',
    '13': 'ccRCC',
    '14': 'Macrophage',
    '15': 'Henle',
    '16': 'Myofibroblast',
    '17': 'NK',
    '18': 'glomerular_endothelial',
    '19': 'Macrophage',
    '20': 'CD4_T',
    '21': 'Myofibroblast',
    '22': 'CD4_T',
    '23': 'Macrophage',
    '24': 'B_cell',
    '25': 'Mast_cell',
    '26': 'collecting_duct',
    '27': 'Erythroid',
    '28': 'proximal',
    '29': 'Podocyte',
}
adata.obs['anno'] = adata.obs['leiden'].astype('str').map(annotation_map)

# before DE, first get rid of MT gene
is_mt = np.logical_not(adata.raw.var['mt'].values)
adata.raw = adata.raw[:,is_mt].to_adata()

# different gene analysis, automatic genes
sc.tl.rank_genes_groups(adata,groupby='anno')
sc.pl.rank_genes_groups(adata,save='.pdf')
sc.pl.rank_genes_groups_heatmap(adata,n_genes=3,swap_axes=True,save='.pdf')

# different gene analysis, interesting genes
genes = ['NDUFA5','NDUFAB1','NDUFAF3','NDUFS3','UQCRC1','MRPL12','MRPL14','MRPL38','MRPL47',
        'MRPS10','MRPS18A','MRPS24','MRPS36','MRPS5','HLA-DMB','HLA-DOA','HLA-DPB1','HLA-DQA2','HLA-DQB2','HLA-DRA']
sc.pl.heatmap(adata,var_names=genes,groupby='anno',swap_axes=True)
plt.savefig('./figures/21gene_de.pdf',bbox_inches='tight')
plt.close()

# different state analysis, automatic
for cluster in adata.obs['anno'].astype('category').cat.categories:
    adata_c = adata[adata.obs['anno']==cluster,:]
    try:
        sc.tl.rank_genes_groups(adata_c,groupby='status')
    except:
        print("cluster {} is in one status\n".format(cluster))
        continue
    else:
        sc.pl.rank_genes_groups(adata_c)
        plt.savefig('./figures/{0}_elbow.pdf'.format(cluster),bbox_inches='tight')
        plt.close()
        sc.pl.rank_genes_groups_heatmap(adata_c,n_genes=15,swap_axes=True)
        plt.savefig('./figures/{0}_heatmap.pdf'.format(cluster),bbox_inches='tight')
        plt.close()

# differential state analysis, interesting genes
for cluster in adata.obs['anno'].astype('category').cat.categories:
    adata_c = adata[adata.obs['anno']==cluster,:]
    try:
        sc.pl.heatmap(adata_c,var_names=genes,groupby='status',swap_axes=True)
    except:
        print("cluster {} is in one status\n".format(cluster))
        continue
    else:
        plt.savefig('./figures/{}_21gene_heatmap.pdf'.format(cluster),bbox_inches='tight')
        plt.close()

# only need ccRCC and 23 genes, and clustering on that
adata_i = adata.raw.to_adata()[adata.obs['anno']=='ccRCC',genes]
df = pd.DataFrame(data=adata_i.X.copy().toarray().T,index=adata_i.var_names)
sns.clustermap(data=df,metric='euclidean',cmap='viridis',row_cluster=False)
plt.savefig('./figures/cluserheatmap.pdf',bbox_inches='tight')
plt.close()


#adata.raw.to_adata().write('./RCC2_to_cellxgene.h5ad')

