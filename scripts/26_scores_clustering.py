import pandas as pd
from scipy.cluster import hierarchy
from matplotlib import pyplot as plt
import seaborn as sns

#clustering patients by the value of their disease
#severity scores (like edss, msss, ecc..)

#reading and cleaning data
df = pd.read_excel('./../database/db.xlsx', sheet_name = 'SM NM', usecols = 'A,CB,DC,FG,JG').dropna()
df.columns = ['patient_id', 'classification', 'edss', 'msss', 'brems']

#selecting RR patients
df = df[df.classification == 3]

#dropping unneeded column and setting appropriate index

df.drop('classification', axis = 1, inplace = True)
df.index = df['patient_id']
df.drop('patient_id', axis = 1, inplace = True)

#clustering dataframe

cluster_col = hierarchy.linkage(df.T, method="ward", metric="euclidean")
cluster_row = hierarchy.linkage(df, method="ward", metric="euclidean")
clusterfig = sns.clustermap(df, row_linkage = cluster_row, col_linkage = cluster_col, cmap = 'mako', yticklabels = True, figsize = (10, len(df)/4))
index_col = clusterfig.dendrogram_col.reordered_ind
index_row = clusterfig.dendrogram_row.reordered_ind
plt.title('RR Clustering')
plt.savefig('../plots/scores_clusters/rr_cluster_scores.png', dpi = 300)
plt.clf()

cluster_col = hierarchy.linkage(df.T, method="ward", metric="euclidean")
cluster_row = hierarchy.linkage(df, method="ward", metric="euclidean")
clusterfig = sns.clustermap(df, row_linkage = cluster_row, col_linkage = cluster_col, cmap = 'mako')
index_col = clusterfig.dendrogram_col.reordered_ind
index_row = clusterfig.dendrogram_row.reordered_ind
plt.title('RR Clustering')
plt.savefig('../plots/scores_clusters/rr_per_slide_scores.png', dpi = 300)
plt.clf()