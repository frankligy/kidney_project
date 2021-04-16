import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np


adata = sc.read('/Volumes/salomonis2/LabFiles/Frank-Li/kidney/RCC1/RCC1_to_cellxgene.h5ad')
bulk = pd.read_csv('/Users/ligk2e/Desktop/tmp/kidney/xena/TCGA.KIRC.sampleMap_HiSeqV2',sep='\t',index_col=0)

sc.pp.highly_variable_genes(adata,flavor='cell_ranger',n_top_genes=1000)
adata = adata[:,adata.var['highly_variable']]


sc_genes = adata.var_names.tolist()
bulk_genes = bulk.index.tolist()
common_genes = list(set(sc_genes).intersection(set(bulk_genes)))


adata = adata[:,common_genes]
bulk = bulk.loc[common_genes,:]

sc_reference = 10**adata.X.toarray() - 1
sc_reference = pd.DataFrame(data=sc_reference,index=adata.obs['anno'].values,columns=adata.var_names)
sc_signature = []  # store multiple series, each of them is a centroid
for cluster in sc_reference.index.unique():
    df = sc_reference.loc[cluster,:]
    sc_signature.append(df.apply(func=np.mean,axis=0))
sc_signature = pd.concat(sc_signature,axis=1)
sc_signature.columns = sc_reference.index.unique()

bulk = 2**bulk - 1

sc_signature.to_csv('/Users/ligk2e/Desktop/tmp/kidney/DeconRNASeq/sc_signature.txt',sep='\t')
bulk.to_csv('/Users/ligk2e/Desktop/tmp/kidney/DeconRNASeq/bulk.txt',sep='\t')


# survival analysis

import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines import CoxPHFitter

final = pd.read_csv('/Users/ligk2e/Desktop/tmp/kidney/DeconRNASeq/final.txt',sep='\t',index_col=0)
survival = pd.read_csv('/Users/ligk2e/Desktop/tmp/kidney/xena/survival.txt',sep='\t',index_col=0)
final.rename(mapper=lambda x:x.replace('.','-'),axis=0,inplace=True)


combine = pd.concat([final,survival],join='inner',axis=1)

sorted_combine = combine.sort_values(by='Tregs')
low = sorted_combine.iloc[:200,:].loc[:,['OS','OS.time']]
high = sorted_combine.iloc[-200:,:].loc[:,['OS','OS.time']]
fig,ax = plt.subplots()
ax.set_ylim(0,1)
for df in [low,high]:
    kmf = KaplanMeierFitter()
    kmf.fit(df['OS.time'],df['OS'])
    kmf.plot_survival_function(ax=ax,ci_show=False,at_risk_counts=False)
current_handles,current_labels = ax.get_legend_handles_labels()
new_labels = ['low_Tregs','high_Tregs']
ax.legend(current_handles,new_labels)
results = logrank_test(low['OS.time'],high['OS.time'],low['OS'],high['OS'])
ax.text(x=1000,y=0.05,s='Log-rank test: p-value is {:.2f}'.format(results.p_value),weight='bold')
plt.savefig('/Users/ligk2e/Desktop/tmp/kidney/xena/Tregs_km.pdf',bbox_inches='tight')

cph = CoxPHFitter()
cph.fit(combine,'OS.time','OS',formula='Macrophage + ccRCC + Tregs + CD8_T + NK')
cph.print_summary()
cph.plot()
plt.savefig('/Users/ligk2e/Desktop/tmp/kidney/xena/cph.pdf',bbox_inches='tight')







