import pandas as pd
from scipy import stats

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

ward_groups = getGroups('./data/cluster_groups/ordered_rr_ward.txt', ward_stops)

#exporting groups in a file

with open('./data/cluster_groups/ward_groups.txt', 'w') as f:
	for group in ward_groups:
		f.write(','.join(str(x) for x in group) + '\n')


group1 = ward_groups[0]
group2 = ward_groups[1]
group3 = ward_groups[2]
group4 = ward_groups[3]
group5 = ward_groups[4]
group6 = ward_groups[5]
group7 = ward_groups[6]

#defining function to get contingency table:
def getContTab(s1, s2):
	succ1 = sum(s1)
	fail1 = len(s1) - succ1
	succ2 = sum(s2)
	fail2 = len(s2) - succ2
	tab = [[succ1, succ2], [fail1, fail2]]
	return tab

def writeStat(matrix, filename):
	with open('./data/clusterbin/{}.tsv'.format(filename), 'w') as f:
		f.write('\tg1\tg2\tg3\tg4\tg5\tg6\tg7\n')
		for i, row in enumerate(matrix, start = 1):
			f.write('g{}\t'.format(i) + '\t'.join(row) + '\n')


####################################### Sesso #######################################

df = pd.read_excel('./../database/db.xlsx', sheet_name = 'SM NM', usecols = 'A,CK')
df.columns = ['patient_id', 'gender']
df1 = df[df['patient_id'].isin(group1)].dropna()
df2 = df[df['patient_id'].isin(group2)].dropna()
df3 = df[df['patient_id'].isin(group3)].dropna()
df4 = df[df['patient_id'].isin(group4)].dropna()
df5 = df[df['patient_id'].isin(group5)].dropna()
df6 = df[df['patient_id'].isin(group6)].dropna()
df7 = df[df['patient_id'].isin(group7)].dropna()
s1 = list(df1['gender'].astype(int))
s2 = list(df2['gender'].astype(int))
s3 = list(df3['gender'].astype(int))
s4 = list(df4['gender'].astype(int))
s5 = list(df5['gender'].astype(int))
s6 = list(df6['gender'].astype(int))
s7 = list(df7['gender'].astype(int))
liste_a = [s1, s2, s3, s4, s5, s6, s7]
liste_b = [s1, s2, s3, s4, s5, s6, s7]
matrix = []

for lista1 in liste_a:
	pvals = []
	for lista2 in liste_b:
		table = getContTab(lista1, lista2)
		pvals.append(str(stats.fisher_exact(table)[1]))
	matrix.append(pvals)
writeStat(matrix, 'sex')

####################################### Nuove Lesioni #######################################

df = pd.read_excel('./../database/db.xlsx', sheet_name = 'SM NM', usecols = 'A,HY')
df.columns = ['patient_id', 'new_les']
df1 = df[df['patient_id'].isin(group1)].dropna()
df2 = df[df['patient_id'].isin(group2)].dropna()
df3 = df[df['patient_id'].isin(group3)].dropna()
df4 = df[df['patient_id'].isin(group4)].dropna()
df5 = df[df['patient_id'].isin(group5)].dropna()
df6 = df[df['patient_id'].isin(group6)].dropna()
df7 = df[df['patient_id'].isin(group7)].dropna()
s1 = list(df1['new_les'].astype(int))
s2 = list(df2['new_les'].astype(int))
s3 = list(df3['new_les'].astype(int))
s4 = list(df4['new_les'].astype(int))
s5 = list(df5['new_les'].astype(int))
s6 = list(df6['new_les'].astype(int))
s7 = list(df7['new_les'].astype(int))
liste_a = [s1, s2, s3, s4, s5, s6, s7]
liste_b = [s1, s2, s3, s4, s5, s6, s7]
matrix = []

for lista1 in liste_a:
	pvals = []
	for lista2 in liste_b:
		table = getContTab(lista1, lista2)
		pvals.append(str(stats.fisher_exact(table)[1]))
	matrix.append(pvals)
writeStat(matrix, 'new_les')

####################################### Gadolinium prima RM #######################################

df = pd.read_excel('./../database/db.xlsx', sheet_name = 'SM NM', usecols = 'A,HZ')
df.columns = ['patient_id', 'gad_mr1']
df1 = df[df['patient_id'].isin(group1)].dropna()
df2 = df[df['patient_id'].isin(group2)].dropna()
df3 = df[df['patient_id'].isin(group3)].dropna()
df4 = df[df['patient_id'].isin(group4)].dropna()
df5 = df[df['patient_id'].isin(group5)].dropna()
df6 = df[df['patient_id'].isin(group6)].dropna()
df7 = df[df['patient_id'].isin(group7)].dropna()
s1 = list(df1['gad_mr1'].astype(int))
s2 = list(df2['gad_mr1'].astype(int))
s3 = list(df3['gad_mr1'].astype(int))
s4 = list(df4['gad_mr1'].astype(int))
s5 = list(df5['gad_mr1'].astype(int))
s6 = list(df6['gad_mr1'].astype(int))
s7 = list(df7['gad_mr1'].astype(int))
liste_a = [s1, s2, s3, s4, s5, s6, s7]
liste_b = [s1, s2, s3, s4, s5, s6, s7]
matrix = []

for lista1 in liste_a:
	pvals = []
	for lista2 in liste_b:
		table = getContTab(lista1, lista2)
		pvals.append(str(stats.fisher_exact(table)[1]))
	matrix.append(pvals)
writeStat(matrix, 'gad_mr1')

####################################### Bande #######################################

df = pd.read_excel('./../database/db.xlsx', sheet_name = 'SM NM', usecols = 'A,DF')
df.columns = ['patient_id', 'bands']
df1 = df[df['patient_id'].isin(group1)].dropna()
df2 = df[df['patient_id'].isin(group2)].dropna()
df3 = df[df['patient_id'].isin(group3)].dropna()
df4 = df[df['patient_id'].isin(group4)].dropna()
df5 = df[df['patient_id'].isin(group5)].dropna()
df6 = df[df['patient_id'].isin(group6)].dropna()
df7 = df[df['patient_id'].isin(group7)].dropna()
s1 = list(df1['bands'].astype(int))
s2 = list(df2['bands'].astype(int))
s3 = list(df3['bands'].astype(int))
s4 = list(df4['bands'].astype(int))
s5 = list(df5['bands'].astype(int))
s6 = list(df6['bands'].astype(int))
s7 = list(df7['bands'].astype(int))
liste_a = [s1, s2, s3, s4, s5, s6, s7]
liste_b = [s1, s2, s3, s4, s5, s6, s7]
matrix = []

for lista1 in liste_a:
	pvals = []
	for lista2 in liste_b:
		table = getContTab(lista1, lista2)
		pvals.append(str(stats.fisher_exact(table)[1]))
	matrix.append(pvals)
writeStat(matrix, 'bands')