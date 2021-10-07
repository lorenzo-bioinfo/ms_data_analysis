import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

#defining function to get groups

def getGroups(file, stops):
	ordered_ids = []
	with open(file) as f:
		for line in f:
			ordered_ids.append(int(line.strip()))
	print(ordered_ids)
	mask = []
	for i in range(len(ordered_ids)):
		if ordered_ids[i] in stops:
			mask.append(-1)
		else:
			mask.append(0)
	print(mask)
	groups = []
	gp = []
	for i in range(len(mask)):
		if mask[i] == 0:
			gp.append(ordered_ids[i])
		else:
			groups.append(gp)
			gp = [ordered_ids[i]]
	groups.append(gp)
	return groups


#extracting groups from cytokines based clustering (ward linkage method)

ward_stops = [1672, 1474, 484, 1577, 1088, 698]

ward_groups = getGroups('../data/cluster_groups/ordered_rr_ward.txt', ward_stops)

group1 = ward_groups[0] #red
group2 = ward_groups[1] #blue
group3 = ward_groups[2] #green
group4 = ward_groups[3] #purple
group5 = ward_groups[4] #orange
group6 = ward_groups[5] #black
group7 = ward_groups[6] #brown

df = pd.read_excel('../../database/db.xlsx', sheet_name = 'SM NM', usecols = 'A,CJ')
df.columns = ['patient_id', 'vitd']
df1 = df[df['patient_id'].isin(group1)].dropna()
df2 = df[df['patient_id'].isin(group2)].dropna()
df3 = df[df['patient_id'].isin(group3)].dropna()
df4 = df[df['patient_id'].isin(group4)].dropna()
df5 = df[df['patient_id'].isin(group5)].dropna()
df6 = df[df['patient_id'].isin(group6)].dropna()
df7 = df[df['patient_id'].isin(group7)].dropna()
s1 = list(df1['vitd'].astype(float))
s2 = list(df2['vitd'].astype(float))
s3 = list(df3['vitd'].astype(float))
s4 = list(df4['vitd'].astype(float))
s5 = list(df5['vitd'].astype(float))
s6 = list(df6['vitd'].astype(float))
s7 = list(df7['vitd'].astype(float))

vectors = [s1, s2, s3, s4, s5, s6, s7]
colors = sns.color_palette('Set1')[0:6]
for i in range(len(colors)):
	sns.kdeplot(vectors[i], color = colors[i])
plt.legend('G1,G2,G3,G4,G5,G6,G7'.split(','))
plt.xlabel('Vitamin D3 Levels')
plt.title('Quantiles Clustering - VitD')
plt.savefig('../../plots/vit_d/quantiles.png', dpi = 300)
plt.clf()

ward_stops = [548, 1394, 1121, 1693, 1107]

ward_groups = getGroups('./../data/cluster_groups/ordered_rr_norm.txt', ward_stops)
group1 = ward_groups[0] #red
group2 = ward_groups[1] #blue
group3 = ward_groups[2] #green
group4 = ward_groups[3] #purple
group5 = ward_groups[4] #orange
group6 = ward_groups[5] #black
df1 = df[df['patient_id'].isin(group1)].dropna()
df2 = df[df['patient_id'].isin(group2)].dropna()
df3 = df[df['patient_id'].isin(group3)].dropna()
df4 = df[df['patient_id'].isin(group4)].dropna()
df5 = df[df['patient_id'].isin(group5)].dropna()
df6 = df[df['patient_id'].isin(group6)].dropna()
s1 = list(df1['vitd'].astype(float))
s2 = list(df2['vitd'].astype(float))
s3 = list(df3['vitd'].astype(float))
s4 = list(df4['vitd'].astype(float))
s5 = list(df5['vitd'].astype(float))
s6 = list(df6['vitd'].astype(float))
vectors = [s1, s2, s3, s4, s5, s6]
colors = sns.color_palette('Set1')[0:5]
for i in range(len(colors)):
	sns.kdeplot(vectors[i], color = colors[i])
plt.legend('G1,G2,G3,G4,G5,G6'.split(','))
plt.xlabel('Vitamin D3 Levels')
plt.title('Normalized Clustering - VitD')
plt.savefig('../../plots/vit_d/normalized.png', dpi = 300)
plt.clf()