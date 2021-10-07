import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
df = pd.read_excel('../database/db.xlsx', sheet_name = 'SM NM', usecols = 'A,F:AF,CB,DD')
df.columns = ('patient_id,IL1B,IL2,IL4,IL5,IL6,IL7,CXCL8,IL10,IL12B,IL13,IL17A,CSF3,CSF2,IFNG,CCL2,CCL4,TNF,IL1RN,IL9,IL15,CCL11,FGF2,CXCL10,PDGFB,CCL5,VEGFA,CCL3,class,cort').split(',')
cyt_list = 'IL1B,IL2,IL4,IL5,IL6,IL7,CXCL8,IL10,IL12B,IL13,IL17A,CSF3,CSF2,IFNG,CCL2,CCL4,TNF,IL1RN,IL9,IL15,CCL11,FGF2,CXCL10,PDGFB,CCL5,VEGFA,CCL3'.split(',')
#dividing df in each class df
df['cort'].fillna(1, inplace = True)
ctrl_df = df[(df['class'] == 6) & (df['cort'] == 1)]
pp_df = df[(df['class'] == 5) & (df['cort'] == 1)]
sp_df = df[(df['class'] == 4) & (df['cort'] == 1)]
rr_df = df[(df['class'] == 3) & (df['cort'] == 1)]

print(ctrl_df, pp_df, sp_df, rr_df)

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

print(len(rr_values), len(pp_values), len(sp_values))
#using the generated arrays to create a df for each class
pp_toclust = pd.DataFrame(pp_values).transpose()
sp_toclust = pd.DataFrame(sp_values).transpose()
rr_toclust = pd.DataFrame(rr_values).transpose()
#assigning cytokine names to corresponding columns
pp_toclust.columns = cyt_list
sp_toclust.columns = cyt_list
rr_toclust.columns = cyt_list
#mapping clusters

sns.clustermap(pp_toclust, method='average', metric='euclidean')
plt.title('PP Clustering')
plt.savefig('../plots/pp_cluster.png', dpi = 300)
plt.clf()

sns.clustermap(sp_toclust, method='average', metric='euclidean')
plt.title('SP Clustering')
plt.savefig('../plots/sp_cluster.png', dpi = 300)
plt.clf()

sns.clustermap(rr_toclust, method='average', metric='euclidean')
plt.title('RR Clustering')
plt.savefig('../plots/rr_cluster.png', dpi = 300)
plt.clf()