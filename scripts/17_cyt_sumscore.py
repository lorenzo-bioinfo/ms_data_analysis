import pandas as pd
from scipy.cluster import hierarchy
from matplotlib import pyplot as plt
import seaborn as sns

cyt_list = 'IL1B,IL2,IL4,IL5,IL6,IL7,CXCL8,IL10,IL12B,IL13,IL17A,CSF3,CSF2,IFNG,CCL2,CCL4,TNF,IL1RN,IL9,IL15,CCL11,FGF2,CXCL10,PDGFB,CCL5,VEGFA,CCL3'.split(',')

#getting df from csv with pre-PL therapy infos

df = pd.read_excel('../database/db.xlsx', sheet_name = 'SM NM', usecols = 'A,F:AF,CB,DD')
df.columns = ('patient_id,IL1B,IL2,IL4,IL5,IL6,IL7,CXCL8,IL10,IL12B,IL13,IL17A,CSF3,CSF2,IFNG,CCL2,CCL4,TNF,IL1RN,IL9,IL15,CCL11,FGF2,CXCL10,PDGFB,CCL5,VEGFA,CCL3,class_int,cort_therapy').split(',')

#This lines were useless, keeping them here fo precaution's sake
#replacing NaN values with 1 (no therapy)
#df['cort_therapy'].fillna(1, inplace = True)

#creating sub-dfs for each class (no therapy before pl)
df_ctrl = df[(df['class_int'] == 6)]
df_ctrl_filled = df_ctrl.fillna({'cort_therapy' : 1})
df_ctrl_noc = df_ctrl_filled[df_ctrl_filled['cort_therapy'] == 1].dropna()
df_pp = df[(df['class_int'] == 5)]
df_pp_filled = df_pp.fillna({'cort_therapy' : 1})
df_pp_noc = df_pp_filled[df_pp_filled['cort_therapy'] == 1].dropna()
df_sp = df[(df['class_int'] == 4)]
df_sp_filled = df_sp.fillna({'cort_therapy' : 1})
df_sp_noc = df_sp_filled[df_sp_filled['cort_therapy'] == 1].dropna()
df_rr = df[(df['class_int'] == 3)]
df_rr_filled = df_rr.fillna({'cort_therapy' : 1})
df_rr_noc = df_rr_filled[df_rr_filled['cort_therapy'] == 1].dropna()
with open('./data/rr_noc_list.txt', 'w') as f:
	lista_rr = list(df_rr_noc['patient_id'].astype(str))
	f.write(','.join(lista_rr))


print('######################################CONTROLLO######################################')
print(df_ctrl_noc) #138
print('######################################PP######################################')
print(df_pp_noc) #44
print('######################################SP######################################')
print(df_sp_noc) #16
print('######################################RR######################################')
print(df_rr_noc) #230

#dropping cortisonic therapy from df as it is useless now
ctrl_df = df_ctrl_noc.drop(columns = 'cort_therapy')
pp_df = df_pp_noc.drop(columns = 'cort_therapy')
sp_df = df_sp_noc.drop(columns = 'cort_therapy')
rr_df = df_rr_noc.drop(columns = 'cort_therapy')

#for each cytokine, defining 10 quantiles.
#These will be used to calculate the cumulative score for
#cytokines classes for each patient
points = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
cyt_quantiles = {}
#based on quartiles
for cyt in cyt_list:
	ctrl_values = ctrl_df[cyt].astype(float)
	quantiles = []
	for point in points:
		quantiles.append(ctrl_values.quantile(point))
	cyt_quantiles[cyt] = quantiles

#defining cytokines groups

cyt_groups = [['CXCL8', 'CCL4', 'CXCL10', 'CCL5'], ['IL12B', 'PDGFB', 'FGF2', 'VEGFA'], ['IFNG', 'IL2', 'IL7', 'IL15', 'CSF3', 'CSF2', 'IL10', 'IL17A', 'IL4', 'IL9', 'IL5', 'IL13'], ['TNF', 'IL6', 'IL1B', 'IL1RN', 'CCL3', 'CCL2', 'CCL11']]

#calculating cumulative scores per cytokine class across RR patients

pat_scores = {}
rr_idx = list(df_rr_noc['patient_id'])

for pat_id in rr_idx:
	pat = df_rr_noc[df_rr_noc['patient_id'] == pat_id]
	pat_scores_list = []
	for group in cyt_groups:
		group_score = []
		for cyt in group:
			added = False
			for quantile in cyt_quantiles[cyt]:
				if list(pat[cyt].astype(float))[0] < quantile:
					group_score.append(points[cyt_quantiles[cyt].index(quantile)] - 0.1)
					added = True
					break
			if not added:
				group_score.append(0.9)
		pat_scores_list.append(sum(group_score) / len(group))
	pat_scores[pat_id] = pat_scores_list

#using dictionary to create df

patscores_df = pd.DataFrame.from_dict(pat_scores, orient = 'index')
patscores_df.columns = ['Cyt1', 'Cyt2', 'Cyt3', 'Cyt4']
print(patscores_df)

#and now clustering

cluster_col = hierarchy.linkage(patscores_df.T, method="ward", metric="euclidean")
cluster_row = hierarchy.linkage(patscores_df, method="ward", metric="euclidean")
clusterfig = sns.clustermap(patscores_df, row_linkage = cluster_row, col_linkage = cluster_col, cmap = 'mako', yticklabels = True, figsize = (10, len(patscores_df)/4))
index_col = clusterfig.dendrogram_col.reordered_ind
index_row = clusterfig.dendrogram_row.reordered_ind
plt.savefig('../plots/rr_sumscore_cluster_ward.png', dpi = 300)
plt.clf()

cluster_col = hierarchy.linkage(patscores_df.T, method="ward", metric="euclidean")
cluster_row = hierarchy.linkage(patscores_df, method="ward", metric="euclidean")
clusterfig = sns.clustermap(patscores_df, row_linkage = cluster_row, col_linkage = cluster_col, cmap = 'mako')
index_col = clusterfig.dendrogram_col.reordered_ind
index_row = clusterfig.dendrogram_row.reordered_ind
plt.savefig('../plots/rr_sumscore_per_slide_ward.png', dpi = 300)
plt.clf()

#getting list of patient ids (not ordered)
ids_list = list(patscores_df.index)
ids_ordered = []
#using the indexes provided by index_row to obtain
#a list of ordered ids (as in the cluster figure)
for i in index_row:
	ids_ordered.append(int(ids_list[i]))

with open('data/cluster_groups/ordered_rr_interactors.txt', 'w') as f:
	for idn in ids_ordered:
		f.write(str(idn) + '\n')