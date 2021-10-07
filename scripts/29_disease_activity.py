import pandas as pd
from scipy import stats

#Trying to see if there is any difference in disease activity
#between groups of patients isolated with different methods

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

#defining functions to perform fisher's exact test

def getContTab(s1, s2):
	succ1 = sum(s1)
	fail1 = len(s1) - succ1
	succ2 = sum(s2)
	fail2 = len(s2) - succ2
	tab = [[succ1, succ2], [fail1, fail2]]
	return tab

def writeStat(matrix, filename):
	with open('./data/activity/{}.tsv'.format(filename), 'w') as f:
		f.write('\tg1\tg2\tg3\tg4\tg5\tg6\tg7\n')
		for i, row in enumerate(matrix, start = 1):
			f.write('g{}\t'.format(i) + '\t'.join(row) + '\n')

#reading data about disease activity:
	#act = activity in general
	#clin = clinical activity
	#rad = radiological activity

df = pd.read_excel('../database/db.xlsx', sheet_name = 'SM NM', usecols = 'A,CN,CO,CP')
df.columns = ['patient_id', 'act', 'clin', 'rad']
df.replace(to_replace = ['NA', 'NO', 2, 3], value = [0, 0, 1, 1], inplace = True)
df = df.dropna()
print(df)

############## Quartiles cytokines clustering ##############

quart_stops = [1672, 1474, 484, 1577, 1088, 698]

quart_groups = getGroups('./data/cluster_groups/ordered_rr_ward.txt', quart_stops)

group1 = quart_groups[0]
group2 = quart_groups[1]
group3 = quart_groups[2]
group4 = quart_groups[3]
group5 = quart_groups[4]
group6 = quart_groups[5]
group7 = quart_groups[6]

df1 = df[df['patient_id'].isin(group1)].dropna()
df2 = df[df['patient_id'].isin(group2)].dropna()
df3 = df[df['patient_id'].isin(group3)].dropna()
df4 = df[df['patient_id'].isin(group4)].dropna()
df5 = df[df['patient_id'].isin(group5)].dropna()
df6 = df[df['patient_id'].isin(group6)].dropna()
df7 = df[df['patient_id'].isin(group7)].dropna()
s1 = list(df1['act'].astype(int))
s2 = list(df2['act'].astype(int))
s3 = list(df3['act'].astype(int))
s4 = list(df4['act'].astype(int))
s5 = list(df5['act'].astype(int))
s6 = list(df6['act'].astype(int))
s7 = list(df7['act'].astype(int))

liste_a = [s1, s2, s3, s4, s5, s6, s7]
liste_b = [s1, s2, s3, s4, s5, s6, s7]
matrix = []
for lista1 in liste_a:
	pvals = []
	for lista2 in liste_b:
		table = getContTab(lista1, lista2)
		pvals.append(str(stats.fisher_exact(table)[1]))
	matrix.append(pvals)
writeStat(matrix, 'quart_act')



s1 = list(df1['clin'].astype(int))
s2 = list(df2['clin'].astype(int))
s3 = list(df3['clin'].astype(int))
s4 = list(df4['clin'].astype(int))
s5 = list(df5['clin'].astype(int))
s6 = list(df6['clin'].astype(int))
s7 = list(df7['clin'].astype(int))

liste_a = [s1, s2, s3, s4, s5, s6, s7]
liste_b = [s1, s2, s3, s4, s5, s6, s7]
matrix = []
for lista1 in liste_a:
	pvals = []
	for lista2 in liste_b:
		table = getContTab(lista1, lista2)
		pvals.append(str(stats.fisher_exact(table)[1]))
	matrix.append(pvals)
writeStat(matrix, 'quart_clin')



s1 = list(df1['rad'].astype(int))
s2 = list(df2['rad'].astype(int))
s3 = list(df3['rad'].astype(int))
s4 = list(df4['rad'].astype(int))
s5 = list(df5['rad'].astype(int))
s6 = list(df6['rad'].astype(int))
s7 = list(df7['rad'].astype(int))

liste_a = [s1, s2, s3, s4, s5, s6, s7]
liste_b = [s1, s2, s3, s4, s5, s6, s7]
matrix = []
for lista1 in liste_a:
	pvals = []
	for lista2 in liste_b:
		table = getContTab(lista1, lista2)
		pvals.append(str(stats.fisher_exact(table)[1]))
	matrix.append(pvals)
writeStat(matrix, 'quart_rad')


############## Sumscore cytokines values clustering ##############

sum_stops = [386, 917, 957, 1463]

sum_groups = getGroups('./data/cluster_groups/ordered_rr_interactors.txt', sum_stops)

group1 = sum_groups[0]
group2 = sum_groups[1]
group3 = sum_groups[2]
group4 = sum_groups[3]
group5 = sum_groups[4]

df1 = df[df['patient_id'].isin(group1)].dropna()
df2 = df[df['patient_id'].isin(group2)].dropna()
df3 = df[df['patient_id'].isin(group3)].dropna()
df4 = df[df['patient_id'].isin(group4)].dropna()
df5 = df[df['patient_id'].isin(group5)].dropna()
s1 = list(df1['act'].astype(int))
s2 = list(df2['act'].astype(int))
s3 = list(df3['act'].astype(int))
s4 = list(df4['act'].astype(int))
s5 = list(df5['act'].astype(int))

liste_a = [s1, s2, s3, s4, s5]
liste_b = [s1, s2, s3, s4, s5]
matrix = []
for lista1 in liste_a:
	pvals = []
	for lista2 in liste_b:
		table = getContTab(lista1, lista2)
		pvals.append(str(stats.fisher_exact(table)[1]))
	matrix.append(pvals)
writeStat(matrix, 'sumscore_act')


s1 = list(df1['clin'].astype(int))
s2 = list(df2['clin'].astype(int))
s3 = list(df3['clin'].astype(int))
s4 = list(df4['clin'].astype(int))
s5 = list(df5['clin'].astype(int))

liste_a = [s1, s2, s3, s4, s5]
liste_b = [s1, s2, s3, s4, s5]
matrix = []
for lista1 in liste_a:
	pvals = []
	for lista2 in liste_b:
		table = getContTab(lista1, lista2)
		pvals.append(str(stats.fisher_exact(table)[1]))
	matrix.append(pvals)
writeStat(matrix, 'sumscore_clin')


s1 = list(df1['rad'].astype(int))
s2 = list(df2['rad'].astype(int))
s3 = list(df3['rad'].astype(int))
s4 = list(df4['rad'].astype(int))
s5 = list(df5['rad'].astype(int))

liste_a = [s1, s2, s3, s4, s5]
liste_b = [s1, s2, s3, s4, s5]
matrix = []
for lista1 in liste_a:
	pvals = []
	for lista2 in liste_b:
		table = getContTab(lista1, lista2)
		pvals.append(str(stats.fisher_exact(table)[1]))
	matrix.append(pvals)
writeStat(matrix, 'sumscore_rad')


############## Normalized cytokines values clustering ##############

norm_stops = [548, 1394, 1121, 1693, 1107]

norm_groups = getGroups('./data/cluster_groups/ordered_rr_norm.txt', norm_stops)

group1 = norm_groups[0]
group2 = norm_groups[1]
group3 = norm_groups[2]
group4 = norm_groups[3]
group5 = norm_groups[4]
group6 = norm_groups[5]


df1 = df[df['patient_id'].isin(group1)].dropna()
df2 = df[df['patient_id'].isin(group2)].dropna()
df3 = df[df['patient_id'].isin(group3)].dropna()
df4 = df[df['patient_id'].isin(group4)].dropna()
df5 = df[df['patient_id'].isin(group5)].dropna()
df6 = df[df['patient_id'].isin(group6)].dropna()
s1 = list(df1['act'].astype(int))
s2 = list(df2['act'].astype(int))
s3 = list(df3['act'].astype(int))
s4 = list(df4['act'].astype(int))
s5 = list(df5['act'].astype(int))
s6 = list(df6['act'].astype(int))

liste_a = [s1, s2, s3, s4, s5, s6]
liste_b = [s1, s2, s3, s4, s5, s6]
matrix = []
for lista1 in liste_a:
	pvals = []
	for lista2 in liste_b:
		table = getContTab(lista1, lista2)
		pvals.append(str(stats.fisher_exact(table)[1]))
	matrix.append(pvals)
writeStat(matrix, 'norm_act')


s1 = list(df1['clin'].astype(int))
s2 = list(df2['clin'].astype(int))
s3 = list(df3['clin'].astype(int))
s4 = list(df4['clin'].astype(int))
s5 = list(df5['clin'].astype(int))
s6 = list(df6['clin'].astype(int))

liste_a = [s1, s2, s3, s4, s5, s6]
liste_b = [s1, s2, s3, s4, s5, s6]
matrix = []
for lista1 in liste_a:
	pvals = []
	for lista2 in liste_b:
		table = getContTab(lista1, lista2)
		pvals.append(str(stats.fisher_exact(table)[1]))
	matrix.append(pvals)
writeStat(matrix, 'norm_clin')


s1 = list(df1['rad'].astype(int))
s2 = list(df2['rad'].astype(int))
s3 = list(df3['rad'].astype(int))
s4 = list(df4['rad'].astype(int))
s5 = list(df5['rad'].astype(int))
s6 = list(df6['rad'].astype(int))

liste_a = [s1, s2, s3, s4, s5, s6]
liste_b = [s1, s2, s3, s4, s5, s6]
matrix = []
for lista1 in liste_a:
	pvals = []
	for lista2 in liste_b:
		table = getContTab(lista1, lista2)
		pvals.append(str(stats.fisher_exact(table)[1]))
	matrix.append(pvals)
writeStat(matrix, 'norm_rad')