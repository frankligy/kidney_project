#! /data/salomonis2/LabFiles/Frank-Li/citeseq/scanpy_new_env/bin/python3.6

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys


# # convert to 10x mtx format
# barcodes = pd.read_csv('./aat1699_DataS1/tableOfCounts_colLabels.tsv',sep='\t',index_col=0) 
# barcodes['DropletID'].to_csv('./aat1699_DataS1/data/barcodes.tsv',sep='\t',index=None,header=None)
# genes = pd.read_csv('./aat1699_DataS1/tableOfCounts_rowLabels.tsv',sep='\t',index_col=0)
# genes.loc[:,['EnsemblID','Symbol']].to_csv('./aat1699_DataS1/data/genes.tsv',sep='\t',index=None,header=None)

# # scanpy workflow
# adata_all = sc.read_10x_mtx('./aat1699_DataS1/data/')
# adata_all.var['mt'] = adata_all.var_names.str.startswith('MT-')
# sc.pp.calculate_qc_metrics(adata_all,qc_vars=['mt'],percent_top=None,inplace=True,log1p=False)

# # first, only consider RCC patients and samples
# meta = pd.read_csv('./meta_sample.txt',sep='\t')
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

# adata.obs.to_csv('./inspection_obs.txt',sep='\t')
# adata.var.to_csv('./inspection_var.txt',sep='\t')


# sc.pp.filter_cells(adata,min_genes=200)
# sc.pp.filter_cells(adata,max_genes=10000)
# sc.pp.filter_genes(adata,min_cells=1)
# adata = adata[adata.obs['pct_counts_mt']<20,:]

# sc.pp.normalize_total(adata,target_sum=1e4)
# sc.pp.log1p(adata)

# adata_ori = adata.copy()
# adata_combat = adata.copy()
# adata_bbknn = adata.copy()

# #original
# adata_ori.raw = adata_ori
# sc.pp.highly_variable_genes(adata_ori,flavor='cell_ranger',n_top_genes=3000)
# adata_ori = adata_ori[:,adata_ori.var['highly_variable']]
# sc.pp.regress_out(adata_ori,['total_counts'])
# sc.pp.scale(adata_ori,max_value=10)
# sc.tl.pca(adata_ori)
# sc.pp.neighbors(adata_ori,n_pcs=50,n_neighbors=30)
# sc.tl.leiden(adata_ori,resolution=1)
# sc.tl.umap(adata_ori)
# adata_ori.write('./adata_ori.h5ad')


# # combat
# adata_combat.raw = adata_combat
# sc.pp.highly_variable_genes(adata_combat,flavor='cell_ranger',n_top_genes=3000)
# adata_combat = adata_combat[:,adata_combat.var['highly_variable']]
# sc.pp.combat(adata_combat,key='channel')
# sc.pp.regress_out(adata_combat,['total_counts'])
# sc.pp.scale(adata_combat,max_value=10)
# sc.tl.pca(adata_combat)
# sc.pp.neighbors(adata_combat,n_pcs=50,n_neighbors=30)
# sc.tl.leiden(adata_combat,resolution=1)
# sc.tl.umap(adata_combat)
# sc.pl.umap(adata_combat,color=['leiden','patient','channel'],legend_loc='right margin')
# plt.savefig('./umap_check.pdf',bbox_inches='tight')
# plt.close()
# adata_combat.write('./adata_combat.h5ad')

# # BBKNN
# adata_bbknn.raw = adata_bbknn
# sc.pp.highly_variable_genes(adata_bbknn,flavor='cell_ranger',n_top_genes=3000)
# adata_bbknn = adata_bbknn[:,adata_bbknn.var['highly_variable']]
# sc.pp.regress_out(adata_bbknn,['total_counts'])
# sc.pp.scale(adata_bbknn,max_value=10)
# c
# sc.external.pp.bbknn(adata_bbknn,batch_key='channel')
# sc.tl.leiden(adata_bbknn,resolution=1)
# sc.tl.umap(adata_bbknn)
# adata_bbknn.write('./adata_bbknn.h5ad')

# # combine to a h5ad, sent to cellxgene
# adata_ori.obsm['X_umap_combat'] = adata_combat.obsm['X_umap']
# adata_ori.obsm['X_umap_bbknn'] = adata_bbknn.obsm['X_umap']
# print(adata_ori)
# adata_ori.raw.to_adata().write('./to_cellxgene.h5ad')


# map kidney annotation to cellxgene h5ad, select ccRCC out, run PAGA, generate visuals
adata = sc.read('./adata_ori.h5ad')
anno = pd.read_csv('./kidney-YMUSK64P.csv',skiprows=[0,1])
barcode_to_anno = pd.Series(data=anno['annotation'].values,index=anno['index']).to_dict()
adata.obs['anno'] = adata.obs_names.map(barcode_to_anno).values
#adata.raw.to_adata().write('./add_anno_to_cellxgene.h5ad')



# sc.tl.paga(adata,groups='anno')
# sc.pl.paga(adata,color='anno',node_size_scale=3,fontsize=5)
# plt.savefig('./paga_plot.pdf',bbox_inches='tight')
# plt.close()

# DE in normal_CD8 and tumor_CD8
sc.tl.rank_genes_groups(adata,groupby='anno',groups=['tumor_CD8'],reference='normal_CD8')
sc.pl.rank_genes_groups(adata,groups=['tumor_CD8'],save='.pdf')
sc.pl.rank_genes_groups_violin(adata,groups=['tumor_CD8'],n_genes=8,strip=False,save='.pdf')

# DE in all
del adata.uns['rank_genes_groups']
sc.tl.rank_genes_groups(adata,groupby='anno')
sc.pl.rank_genes_groups_stacked_violin(adata,n_genes=3,save='.pdf')
sc.pl.rank_genes_groups_dotplot(adata,n_genes=3,save='.pdf')
sc.pl.rank_genes_groups_matrixplot(adata,n_genes=3,save='.pdf')
sc.pl.rank_genes_groups_heatmap(adata,n_genes=3,swap_axes=True,save='.pdf')
sc.pl.rank_genes_groups_tracksplot(adata,n_genes=3,save='.png')





'''
Only focus on ccRCC
'''
adata_t = adata[adata.obs['anno']=='ccRCC',:]
sc.tl.pca(adata_t)
sc.pp.neighbors(adata_t,n_pcs=50,n_neighbors=30)
sc.tl.leiden(adata_t,resolution=0.5)
sc.tl.umap(adata_t)
sc.pl.umap(adata_t,color='leiden')
# plt.savefig('./ccRCC_umap.pdf',bbox_inches='tight')
# plt.close()
# adata_t.raw.to_adata().write('./ccRCC_to_cellxgene.h5ad')


sc.tl.paga(adata_t,groups='leiden')
sc.pl.paga(adata_t,color='leiden',node_size_scale=4,fontsize=3)
plt.savefig('./paga_plot.pdf',bbox_inches='tight')
plt.close()

sc.tl.diffmap(adata_t)
adata_t.uns['iroot'] = np.flatnonzero(adata_t.obs['leiden'] == '0')[0]
sc.tl.dpt(adata_t)
sc.pl.umap(adata_t,color=['dpt_pseudotime'])
plt.savefig('./umap_embed_dpt.pdf',bbox_inches='tight')
plt.close()

# draw paga_path
gene_names = ['CA9','NDUFA4L2','HIF1A','VHL','MTOR','HIF3A','PBRM1','SETD2','BAP1']
path = ['0','5','3','4','6','2','1']
adata_t = adata_t.raw.to_adata()
print(adata_t)
print(adata_t.obs['leiden'].cat.categories)
sc.pl.paga_path(adata_t,nodes=path,keys=gene_names,use_raw=False)
plt.savefig('./paga_path.pdf',bbox_inches='tight')
plt.close()











