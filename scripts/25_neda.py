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
	with open('./data/neda_stats/{}.tsv'.format(filename), 'w') as f:
		f.write('\tg1\tg2\tg3\tg4\tg5\tg6\tg7\n')
		for i, row in enumerate(matrix, start = 1):
			f.write('g{}\t'.format(i) + '\t'.join(row) + '\n')

#reading and cleaning data

df = pd.read_excel('./../database/db.xlsx', sheet_name = 'SM NM', usecols = 'A,JC').dropna()
df.columns = ['patient_id', 'neda']
df = df[df.neda != '?']
df = df[df.neda != 'x']
df = df[df.neda != 0.509714]

#testing differences between quartiles groups

df1 = df[df['patient_id'].isin(group1)].dropna()
df2 = df[df['patient_id'].isin(group2)].dropna()
df3 = df[df['patient_id'].isin(group3)].dropna()
df4 = df[df['patient_id'].isin(group4)].dropna()
df5 = df[df['patient_id'].isin(group5)].dropna()
df6 = df[df['patient_id'].isin(group6)].dropna()
df7 = df[df['patient_id'].isin(group7)].dropna()
s1 = list(df1['neda'].astype(int))
s2 = list(df2['neda'].astype(int))
s3 = list(df3['neda'].astype(int))
s4 = list(df4['neda'].astype(int))
s5 = list(df5['neda'].astype(int))
s6 = list(df6['neda'].astype(int))
s7 = list(df7['neda'].astype(int))
liste_a = [s1, s2, s3, s4, s5, s6, s7]
liste_b = [s1, s2, s3, s4, s5, s6, s7]
matrix = []

for lista1 in liste_a:
	pvals = []
	for lista2 in liste_b:
		table = getContTab(lista1, lista2)
		pvals.append(str(stats.fisher_exact(table)[1]))
	matrix.append(pvals)
writeStat(matrix, 'quart')

#testing differences between norm groups

def writeStat(matrix, filename):
	with open('./data/neda_stats/{}.tsv'.format(filename), 'w') as f:
		f.write('\tg1\tg2\tg3\tg4\tg5\tg6\n')
		for i, row in enumerate(matrix, start = 1):
			f.write('g{}\t'.format(i) + '\t'.join(row) + '\n')

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

s1 = list(df1['neda'].astype(int))
s2 = list(df2['neda'].astype(int))
s3 = list(df3['neda'].astype(int))
s4 = list(df4['neda'].astype(int))
s5 = list(df5['neda'].astype(int))
s6 = list(df6['neda'].astype(int))

liste_a = [s1, s2, s3, s4, s5, s6]
liste_b = [s1, s2, s3, s4, s5, s6]
matrix = []

for lista1 in liste_a:
	pvals = []
	for lista2 in liste_b:
		table = getContTab(lista1, lista2)
		pvals.append(str(stats.fisher_exact(table)[1]))
	matrix.append(pvals)
writeStat(matrix, 'norm')

#testing differences between sumscore groups

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

s1 = list(df1['neda'].astype(int))
s2 = list(df2['neda'].astype(int))
s3 = list(df3['neda'].astype(int))
s4 = list(df4['neda'].astype(int))
s5 = list(df5['neda'].astype(int))

liste_a = [s1, s2, s3, s4, s5]
liste_b = [s1, s2, s3, s4, s5]
matrix = []

def writeStat(matrix, filename):
	with open('./data/neda_stats/{}.tsv'.format(filename), 'w') as f:
		f.write('\tg1\tg2\tg3\tg4\tg5\n')
		for i, row in enumerate(matrix, start = 1):
			f.write('g{}\t'.format(i) + '\t'.join(row) + '\n')

for lista1 in liste_a:
	pvals = []
	for lista2 in liste_b:
		table = getContTab(lista1, lista2)
		pvals.append(str(stats.fisher_exact(table)[1]))
	matrix.append(pvals)
writeStat(matrix, 'sumscore')