import pandas as pd
from sklearn.preprocessing import MinMaxScaler
from scipy.cluster import hierarchy
from matplotlib import pyplot as plt
import seaborn as sns

#clustering patients by the value of their disease
#severity scores (like edss, msss, ecc..) after
#normalizing the values

#reading and cleaning data
df = pd.read_excel('./../database/db.xlsx', sheet_name = 'SM NM', usecols = 'A,CB,DC,FG,JG').dropna()
df.columns = ['patient_id', 'classification', 'edss', 'msss', 'brems']

#selecting RR patients
df = df[df.classification == 3]

#dropping unneeded column and setting appropriate index

df.drop('classification', axis = 1, inplace = True)
df.index = df['patient_id']
df.drop('patient_id', axis = 1, inplace = True)

#applying MinMaxScaler to normalize scores data

pat_id = list(df.index)
scores_list = ['edss', 'msss', 'brems']
scores_matrix = []
for score in scores_list:
	arr = list(df[score].astype(float))
	scores_matrix.append(arr)
df_tonorm = pd.DataFrame(scores_matrix).transpose()
norm = MinMaxScaler().fit(df_tonorm)
df_norm = norm.transform(df_tonorm)
df_norm = pd.DataFrame(df_norm)
df_norm.index = pat_id
df_norm.columns = scores_list

#clustering normalized dataframe

cluster_col = hierarchy.linkage(df_norm.T, method="ward", metric="euclidean")
cluster_row = hierarchy.linkage(df_norm, method="ward", metric="euclidean")
clusterfig = sns.clustermap(df_norm, row_linkage = cluster_row, col_linkage = cluster_col, cmap = 'mako', yticklabels = True, figsize = (10, len(df_norm)/4))
index_col = clusterfig.dendrogram_col.reordered_ind
index_row = clusterfig.dendrogram_row.reordered_ind
plt.title('RR Clustering')
plt.savefig('../plots/scores_clusters/rr_norm_cluster_scores.png', dpi = 300)
plt.clf()

cluster_col = hierarchy.linkage(df_norm.T, method="ward", metric="euclidean")
cluster_row = hierarchy.linkage(df_norm, method="ward", metric="euclidean")
clusterfig = sns.clustermap(df_norm, row_linkage = cluster_row, col_linkage = cluster_col, cmap = 'mako')
index_col = clusterfig.dendrogram_col.reordered_ind
index_row = clusterfig.dendrogram_row.reordered_ind
plt.title('RR Clustering')
plt.savefig('../plots/scores_clusters/rr_norm_per_slide_scores.png', dpi = 300)
plt.clf()

#getting list of patient ids (not ordered)
ids_list = list(df_norm.index)
ids_ordered = []
#using the indexes provided by index_row to obtain
#a list of ordered ids (as in the cluster figure)
for i in index_row:
	ids_ordered.append(ids_list[i])

#exporting ordered patients ids

with open('data/cluster_groups/ordered_rr_sevscores_norm.txt', 'w') as f:
	for idn in ids_ordered:
		f.write(str(idn) + '\n')