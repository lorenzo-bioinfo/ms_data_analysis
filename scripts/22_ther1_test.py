#In this script I will try to assess any possible difference
#in treatment outcome across patients of different groups
#(isolated with hierarchical clustering on quartiles cyt values)

import pandas as pd
from scipy import stats

#importing table with treatment infos

df_rr = pd.read_csv('./checks/delta_pl_ther.tsv', sep = '\t')

#dropping patient where difference between PL and
#first day of first treatment is higher than 90 days
df = df_rr[(abs(df_rr['delta_pl']) < 91)]
print(df)

#selecting patients groups

def getGroups(file, stops):
	ordered_ids = []
	with open(file) as f:
		for line in f:
			ordered_ids.append(int(line.strip()))
	print('NUMERO PAZIENTI: ', len(ordered_ids))
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

int_stops = [1672, 1474, 484, 1577, 1088, 698]

groups = getGroups('./data/cluster_groups/ordered_rr_ward.txt', int_stops)

#creating a different df for each group

group1 = groups[0]
group2 = groups[1]
group3 = groups[2]
group4 = groups[3]
group5 = groups[4]
group6 = groups[5]
group7 = groups[6]

df1 = df[df['patient_id'].isin(group1)]
df2 = df[df['patient_id'].isin(group2)]
df3 = df[df['patient_id'].isin(group3)]
df4 = df[df['patient_id'].isin(group4)]
df5 = df[df['patient_id'].isin(group5)]
df6 = df[df['patient_id'].isin(group6)]
df7 = df[df['patient_id'].isin(group7)]

#defining function to get contingency table:
def getContTab(s1, s2):
	succ1 = sum(s1)
	fail1 = len(s1) - succ1
	succ2 = sum(s2)
	fail2 = len(s2) - succ2
	tab = [[succ1, succ2], [fail1, fail2]]
	return tab

def writeStat(matrix, filename):
	with open('./data/thertest/{}.tsv'.format(filename), 'w') as f:
		f.write('\tg1\tg2\tg3\tg4\tg5\tg6\tg7\n')
		for i, row in enumerate(matrix, start = 1):
			f.write('g{}\t'.format(i) + '\t'.join(row) + '\n')

#testing if patients of different groups have a different response
#to therapy (without distinction of FL vs SL)

def getInterruptionCause(df):
	series = df['re_end1'].fillna('A')
	series.replace(to_replace = 'P', value = 1, inplace = True)
	series.replace(to_replace = 'A', value = 0, inplace = True)
	series.dropna(inplace = True)
	return(list(series))

s1 = getInterruptionCause(df1)
s2 = getInterruptionCause(df2)
s3 = getInterruptionCause(df3)
s4 = getInterruptionCause(df4)
s5 = getInterruptionCause(df5)
s6 = getInterruptionCause(df6)
s7 = getInterruptionCause(df7)

liste_a = [s1, s2, s3, s4, s5, s6, s7]
liste_b = [s1, s2, s3, s4, s5, s6, s7]
matrix = []

for lista1 in liste_a:
	pvals = []
	for lista2 in liste_b:
		table = getContTab(lista1, lista2)
		print(table)
		pvals.append(str(stats.fisher_exact(table)[1]))
	matrix.append(pvals)
writeStat(matrix, 'ther1_end')

#testing if patients of different groups have been assigned
#to different therapy classes (first or second line)

def getTherapyLine(df):
	series = df['TL'].dropna()
	series.replace(to_replace = 'FL', value = 1, inplace = True)
	series.replace(to_replace = 'SL', value = 0, inplace = True)
	series.dropna(inplace = True)
	return(list(series))

s1 = getTherapyLine(df1)
s2 = getTherapyLine(df2)
s3 = getTherapyLine(df3)
s4 = getTherapyLine(df4)
s5 = getTherapyLine(df5)
s6 = getTherapyLine(df6)
s7 = getTherapyLine(df7)

liste_a = [s1, s2, s3, s4, s5, s6, s7]
liste_b = [s1, s2, s3, s4, s5, s6, s7]
matrix = []

for lista1 in liste_a:
	pvals = []
	for lista2 in liste_b:
		table = getContTab(lista1, lista2)
		print(table)
		pvals.append(str(stats.fisher_exact(table)[1]))
	matrix.append(pvals)
writeStat(matrix, 'ther_line')