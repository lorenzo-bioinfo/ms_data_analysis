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

#defining thresholds for low, normal and high level for each cytokine
thres = []
#based on quartiles
for cyt in cyt_list:
	thres.append((ctrl_df[cyt].astype(float).quantile(0.25), ctrl_df[cyt].astype(float).quantile(0.75)))
#associating the threshold values to each cytokine using a dictionary
cyt_thres = dict(zip(cyt_list, thres))
for key in cyt_thres:
	print('{} - {}'.format(key, cyt_thres[key]))
#for each class:
	#setting low values to -1
	#normal values to 0
	#high values to 1
pp_values = []
for cyt in cyt_list:
	arr = list(pp_df[cyt].astype(float))
	to_push = []
	for entry in arr:
		if entry < cyt_thres[cyt][0]:
			entry = -1
		elif entry > cyt_thres[cyt][1]:
			entry = 1
		else:
			entry = 0
		to_push.append(entry)
	pp_values.append(to_push)

pp_toclust = pd.DataFrame(pp_values).transpose()
pp_toclust.index = list(pp_df['patient_id'])
pp_toclust.columns = cyt_list

sp_values = []
for cyt in cyt_list:
	arr = list(sp_df[cyt].astype(float))
	to_push = []
	for entry in arr:
		if entry < cyt_thres[cyt][0]:
			entry = -1
		elif entry > cyt_thres[cyt][1]:
			entry = 1
		else:
			entry = 0
		to_push.append(entry)
	sp_values.append(to_push)

sp_toclust = pd.DataFrame(sp_values).transpose()
sp_toclust.index = list(sp_df['patient_id'])
sp_toclust.columns = cyt_list

rr_values = []
for cyt in cyt_list:
	arr = list(rr_df[cyt].astype(float))
	to_push = []
	for entry in arr:
		if entry < cyt_thres[cyt][0]:
			entry = -1
		elif entry > cyt_thres[cyt][1]:
			entry = 1
		else:
			entry = 0
		to_push.append(entry)
	rr_values.append(to_push)

rr_toclust = pd.DataFrame(rr_values).transpose()
rr_toclust.index = list(rr_df['patient_id'])
rr_toclust.columns = cyt_list
print(len(rr_toclust))

##############################################Recursively removing groups##################################################

with open('./data/cluster_groups/ward_groups.txt', 'r') as f:
	groups = []
	for line in f:
		groups.append(line.strip().split(','))

cyt_ord = 'TNF,IL13,IL9,IL1B,IL7,IL4,VEGFA,IL10,CCL11,IL2,IL12B,IL15,PDGFB,IL6,IL17A,CSF2,FGF2,CCL4,CCL3,IFNG,CCL2,CXCL8,CXCL10,IL5,IL1RN,CSF3,CCL5'.split(',')
rr_forced = pd.DataFrame()

for cyt in cyt_ord:
	rr_forced[cyt] = rr_toclust[cyt]
rr_forced.index = rr_toclust.index

clean_groups = []
for group in groups:
	a = [int(x) for x in group]
	clean_groups.append(a)

clean_groups.pop()

i = 1
for group in clean_groups:
	rr_forced.drop(group, axis = 0, inplace = True)

	cluster_row = hierarchy.linkage(rr_forced, method="ward", metric="euclidean")
	cluster_col = hierarchy.linkage(rr_forced.T, method="ward", metric="euclidean")
	clusterfig = sns.clustermap(rr_forced, row_linkage = cluster_row, col_linkage = cluster_col, yticklabels = True, figsize = (10, len(rr_forced)/4))
	index_row = clusterfig.dendrogram_row.reordered_ind
	plt.title('RR Clustering')
	plt.savefig('../plots/remove_clusters/rr_unforced_ward_noG{}.png'.format(i), dpi = 300)
	plt.clf()

	cluster_row = hierarchy.linkage(rr_forced, method="ward", metric="euclidean")
	cluster_col = hierarchy.linkage(rr_forced.T, method="ward", metric="euclidean")
	clusterfig = sns.clustermap(rr_forced, row_linkage = cluster_row, col_linkage = cluster_col)
	index_row = clusterfig.dendrogram_row.reordered_ind
	plt.title('RR Clustering')
	plt.savefig('../plots/remove_clusters/small_rr_unforced_ward_noG{}.png'.format(i), dpi = 300)
	plt.clf()

	i +=1
