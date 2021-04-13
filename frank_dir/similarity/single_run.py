#! /data/salomonis2/LabFiles/Frank-Li/citeseq/scanpy_new_env/bin/python3.6

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import scrublet as scr
# import mnnpy



# read RCC1 and RCC2
adata1 = sc.read('../RCC1/RCC1_after_umap.h5ad')
adata2 = sc.read('../RCC2/RCC2_after_umap.h5ad')

# put raw back to adata
adata1 = adata1.raw.to_adata()
adata2 = adata2.raw.to_adata()

# get common genes
rcc1_genes = set(adata1.var_names.to_list())
rcc2_genes = set(adata2.var_names.to_list())
common = list(rcc1_genes.intersection(rcc2_genes))

adata1 = adata1[:,common]
adata2 = adata2[:,common]

# # perform batch-effect removal 
# sc.pp.highly_variable_genes(adata1,flavor='cell_ranger',n_top_genes=3000)
# sc.pp.highly_variable_genes(adata2,flavor='cell_ranger',n_top_genes=3000)
# hvg_adata1 = adata1[:,adata1.var['highly_variable']].var_names.to_list()
# hvg_adata2 = adata2[:,adata2.var['highly_variable']].var_names.to_list()
# hvg = list(set(hvg_adata1).intersection(set(hvg_adata2)))


# corrected = mnnpy.mnn_correct(adata1,adata2,var_subset=hvg,do_concatenate=False)
# adata1_post,adata2_post = corrected[0]
# adata1_post.write('./adata1_post_mnn.h5ad')
# adata2_post.write('./adata2_post_mnn.h5ad')

adata1_post = sc.read('./adata1_post_mnn.h5ad')
adata2_post = sc.read('./adata2_post_mnn.h5ad')

# adata12_post = adata1_post.concatenate(adata2_post)
# sc.pp.highly_variable_genes(adata12_post,flavor='cell_ranger',n_top_genes=3000)
# adata12_post = adata12_post[:,adata12_post.var['highly_variable']]
# sc.pp.regress_out(adata12_post,['total_counts'])
# sc.pp.scale(adata12_post,max_value=10)
# sc.tl.pca(adata12_post)
# sc.pp.neighbors(adata12_post)
# sc.tl.umap(adata12_post)
# adata12_post.write('./adata12_post_mnn.h5ad')






# add annotations
annotation_map1 ={
    '0': 'Macrophage',
    '1': 'glomerular_endothelial',
    '2': 'proximal',
    '3': 'glomerular_endothelial',
    '4': 'proximal',
    '5': 'ccRCC',
    '6': 'CD8_T',
    '7': 'NK',
    '8': 'CD4_T',
    '9': 'ccRCC',
    '10': 'CD4_T',
    '11': 'Macrophage',
    '12': 'Henle',
    '13': 'glomerular_endothelial',
    '14': 'Macrophage',
    '15': 'Tregs',
    '16': 'glomerular_endothelial',
    '17': 'collecting_duct',
    '18': 'myofibroblast',
    '19': 'NK',
    '20': 'CD8_T',
    '21': 'podocyte',
    '22': 'Macrophage',
    '23': 'Erythroid',
    '24': 'B_cell',
    '25': 'glomerular_endothelial',
}
adata1.obs['anno'] = adata1.obs['leiden'].astype('str').map(annotation_map1)
adata1_post.obs['anno'] = adata1_post.obs['leiden'].astype('str').map(annotation_map1)


annotation_map2 ={
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
adata2.obs['anno'] = adata2.obs['leiden'].astype('str').map(annotation_map2)
adata2_post.obs['anno'] = adata2_post.obs['leiden'].astype('str').map(annotation_map2)



# only tumor
adata1_c = adata1[adata1.obs['anno']=='ccRCC',:]
adata2_c = adata2[adata2.obs['anno']=='ccRCC',:]
adata1_post_c = adata1_post[adata1_post.obs['anno']=='ccRCC',:]
adata2_post_c = adata2_post[adata2_post.obs['anno']=='ccRCC',:]




# leiden on both
adata1_c.raw = adata1_c
sc.pp.filter_genes(adata1_c,min_cells=1)
sc.pp.highly_variable_genes(adata1_c,flavor='cell_ranger',n_top_genes=3000)
adata1_c = adata1_c[:,adata1_c.var['highly_variable']]
sc.tl.pca(adata1_c)
sc.pp.neighbors(adata1_c)
sc.tl.leiden(adata1_c,resolution=0.5)
sc.tl.umap(adata1_c)
sc.pl.umap(adata1_c,color=['leiden'])
plt.savefig('./RCC1_cancer_umap.pdf',bbox_inches='tight')
plt.close()
sc.tl.rank_genes_groups(adata1_c,groupby='leiden')
sc.pl.rank_genes_groups_heatmap(adata1_c,n_genes=5,swap_axes=True)
plt.savefig('./RCC1_cancer_heatmap.pdf',bbox_inches='tight')
plt.close()

# add the interested genes



adata2_c.raw = adata2_c
sc.pp.filter_genes(adata2_c,min_cells=1)
sc.pp.highly_variable_genes(adata2_c,flavor='cell_ranger',n_top_genes=3000)
adata2_c = adata2_c[:,adata2_c.var['highly_variable']]
sc.tl.pca(adata2_c)
sc.pp.neighbors(adata2_c)
sc.tl.leiden(adata2_c,resolution=0.5)
sc.tl.umap(adata2_c)
sc.pl.umap(adata2_c,color=['leiden'])
plt.savefig('./RCC2_cancer_umap.pdf',bbox_inches='tight')
plt.close()
sc.tl.rank_genes_groups(adata2_c,groupby='leiden')
sc.pl.rank_genes_groups_heatmap(adata2_c,n_genes=5,swap_axes=True)
plt.savefig('./RCC2_cancer_heatmap.pdf',bbox_inches='tight')
plt.close()



# # let's see what are the differentially expressed genes between two tumors
# adata12_c = adata1_c.concatenate(adata2_c)
# sc.tl.rank_genes_groups(adata12_c,groupby='patient')
# sc.pl.rank_genes_groups_heatmap(adata12_c,n_genes=25,swap_axes=True)
# plt.savefig('./DE_between12.pdf',bbox_inches='tight')
# plt.close()
# for patient in adata12_c.obs['patient'].astype('category').cat.categories:
#     df = sc.get.rank_genes_groups_df(adata12_c,group=patient)
#     df.to_csv('de_list_{}.txt'.format(patient),sep='\t',index=None)



# umap ON post-mnn
post = sc.read('./adata12_post_mnn.h5ad')
post1 = post[post.obs['patient']=='RCC1',:]
post1.obs['anno'] = post1.obs['leiden'].astype('str').map(annotation_map1)
post2 = post[post.obs['patient']=='RCC2',:]
post2.obs['anno'] = post2.obs['leiden'].astype('str').map(annotation_map2)

adata1_c.obs.to_csv('check1.txt',sep='\t')

post1_c = post1[post1.obs['anno']=='ccRCC',:]
post1_c.obs['final'] = ['RCC1_'+item for item in adata1_c.obs['leiden']]
post2_c = post2[post2.obs['anno']=='ccRCC',:]
post2_c.obs['final'] = ['RCC2_'+item for item in adata2_c.obs['leiden']]

post12 = post1_c.concatenate(post2_c)
sc.pl.umap(post12,color='final')
plt.savefig('./layout.pdf',bbox_inches='tight')
plt.close()

post12.write('post12_to_cellxgene.h5ad')

# just for batch-effect removal
adata1_post_c.obs['leiden'] = adata1_c.obs['leiden']
adata2_post_c.obs['leiden'] = adata2_c.obs['leiden']
adata1_post_c.raw = adata1_post_c
adata2_post_c.raw = adata2_post_c
adata1_c = adata1_post_c
adata2_c = adata2_post_c





''' train in RCC1, test on RCC2'''
from sklearn.linear_model import LogisticRegression
# train on RCC1
X = adata1_c.raw.X
Y = adata1_c.obs['leiden'].values
model = LogisticRegression(penalty='l1',solver='liblinear',max_iter=100000)
model.fit(X,Y)

# test on RCC2
X_test = adata2_c.raw.X
Y_test = adata2_c.obs['leiden'].values  # (n_cells,)
result = model.predict(X_test)  # (n_cells,)

# construct confusion matrix
mat = np.zeros([len(adata2_c.obs['leiden'].cat.categories),len(adata1_c.obs['leiden'].cat.categories)])  # row is number of categores in RCC2, column is the number in RCC1
for i,cluster in enumerate(adata2_c.obs['leiden'].cat.categories):
    original = Y_test[np.where(Y_test==cluster)[0]]
    prediction = result[np.where(Y_test==cluster)[0]]
    u,c = np.unique(prediction,return_counts=True)
    for j in range(len(u)):
        mat[i,int(u[j])] = int(c[j])/c.sum()
mat_df12 = pd.DataFrame(data=mat,index=['RCC2_'+item for item in adata2_c.obs['leiden'].cat.categories],
                        columns=['RCC1_'+item for item in adata1_c.obs['leiden'].cat.categories])
import seaborn as sns
sns.heatmap(mat_df12,annot=True)
plt.savefig('./RCC1-RCC2.pdf',bbox_inches='tight')
plt.close()

''' train in RCC2, test on RCC1'''
from sklearn.linear_model import LogisticRegression
# train on RCC2
X = adata2_c.raw.X
Y = adata2_c.obs['leiden'].values
model = LogisticRegression(penalty='l1',solver='liblinear',max_iter=100000)
model.fit(X,Y)

# test on RCC1
X_test = adata1_c.raw.X
Y_test = adata1_c.obs['leiden'].values  # (n_cells,)
result = model.predict(X_test)  # (n_cells,)

# construct confusion matrix
mat = np.zeros([len(adata1_c.obs['leiden'].cat.categories),len(adata2_c.obs['leiden'].cat.categories)])  # row is number of categores in RCC2, column is the number in RCC1
for i,cluster in enumerate(adata1_c.obs['leiden'].cat.categories):
    original = Y_test[np.where(Y_test==cluster)[0]]
    prediction = result[np.where(Y_test==cluster)[0]]
    u,c = np.unique(prediction,return_counts=True)
    for j in range(len(u)):
        mat[i,int(u[j])] = int(c[j])/c.sum()
mat_df21 = pd.DataFrame(data=mat,index=['RCC1_'+item for item in adata1_c.obs['leiden'].cat.categories],
                        columns=['RCC2_'+item for item in adata2_c.obs['leiden'].cat.categories])
sns.heatmap(mat_df21,annot=True)
plt.savefig('./RCC2-RCC1.pdf',bbox_inches='tight')
plt.close()




# '''avearge two confusion matrix'''
# mat_df = pd.DataFrame(data=(mat_df21.values + mat_df12.T.values)/2,index=['RCC1_'+item for item in adata1_c.obs['leiden'].cat.categories],
#                         columns=['RCC2_'+item for item in adata2_c.obs['leiden'].cat.categories])
# sns.heatmap(mat_df,annot=True)
# plt.savefig('./RCC.pdf',bbox_inches='tight')
# plt.close()





# '''Some additional gene list'''
# hla_gene = ['HLA-A','HLA-B','HLA-C','HLA-DMA','HLA-DMB','HLA-DOA','HLA-DOB','HLA-DPA1','HLA-DPA2','HLA-DPA3','HLA-DPB1','HLA-DPB2',
#             'HLA-DQA1','HLA-DQA2','HLA-DQB1','HLA-DQB2','HLA-DQB3','HLA-DRA','HLA-DRB1','HLA-DRB2','HLA-DRB3','HLA-DRB4','HLA-DRB5',
#             'HLA-DRB6','HLA-DRB7','HLA-DRB8','HLA-DRB9']
# copper_gene = ['MT-CO1','MT-CO2','MT-CO3']



