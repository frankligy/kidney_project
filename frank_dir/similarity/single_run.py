#! /data/salomonis2/LabFiles/Frank-Li/citeseq/scanpy_new_env/bin/python3.6

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import scrublet as scr
import sys


# read RCC1 and RCC2
adata1 = sc.read('../RCC1/RCC1_after_umap.h5ad')
adata2 = sc.read('../RCC2/RCC2_after_umap.h5ad')

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

# put raw back to adata
adata1 = adata1.raw.to_adata()
adata1 = adata1[:,~adata1.var['mt']]


adata2 = adata2.raw.to_adata()
adata2 = adata2[:,~adata2.var['mt']]


# get common genes and only pull out ccRCC
rcc1_genes = set(adata1.var_names.to_list())
rcc2_genes = set(adata2.var_names.to_list())
common = list(rcc1_genes.intersection(rcc2_genes))

adata1_c = adata1[adata1.obs['anno']=='ccRCC',common]
adata2_c = adata2[adata2.obs['anno']=='ccRCC',common]


# leiden on both
adata1_c.raw = adata1_c
sc.pp.filter_genes(adata1_c,min_cells=1)
sc.pp.highly_variable_genes(adata1_c,flavor='cell_ranger',n_top_genes=3000)
adata1_c = adata1_c[:,adata1_c.var['highly_variable']]
sc.tl.pca(adata1_c)
sc.pp.neighbors(adata1_c)
sc.tl.leiden(adata1_c,resolution=0.5)
sc.tl.umap(adata1_c)
sc.tl.rank_genes_groups(adata1_c,groupby='leiden')
sc.pl.rank_genes_groups_heatmap(adata1_c,n_genes=5,swap_axes=True)
plt.savefig('./RCC1_cancer_heatmap.pdf',bbox_inches='tight')
plt.close()



adata2_c.raw = adata2_c
sc.pp.filter_genes(adata2_c,min_cells=1)
sc.pp.highly_variable_genes(adata2_c,flavor='cell_ranger',n_top_genes=3000)
adata2_c = adata2_c[:,adata2_c.var['highly_variable']]
sc.tl.pca(adata2_c)
sc.pp.neighbors(adata2_c)
sc.tl.leiden(adata2_c,resolution=0.5)
sc.tl.umap(adata2_c)
sc.tl.rank_genes_groups(adata2_c,groupby='leiden')
sc.pl.rank_genes_groups_heatmap(adata2_c,n_genes=5,swap_axes=True)
plt.savefig('./RCC2_cancer_heatmap.pdf',bbox_inches='tight')
plt.close()


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


'''avearge two confusion matrix'''
mat_df = pd.DataFrame(data=(mat_df21.values + mat_df12.T.values)/2,index=['RCC1_'+item for item in adata1_c.obs['leiden'].cat.categories],
                        columns=['RCC2_'+item for item in adata2_c.obs['leiden'].cat.categories])
sns.heatmap(mat_df,annot=True)
plt.savefig('./RCC.pdf',bbox_inches='tight')
plt.close()








