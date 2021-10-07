import pandas as pd
from scipy import stats
from numpy import median

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

for group in ward_groups:
	#print(group)
	print(len(group))

print(ward_groups[0])

#exporting groups in a file

with open('./data/cluster_groups/ward_groups.txt', 'w') as f:
	for group in ward_groups:
		f.write(','.join(str(x) for x in group) + '\n')

#Performing Kolmogorov-Smirnov test on different parameters
#Each group will be tested against all others

def cross_ks(s1, s2, s3, s4, s5, s6, s7, filename):
	print('Testing {}'.format(filename))
	tup1 = (s1, s2, s3, s4, s5, s6, s7)
	tup2 = (s1, s2, s3, s4, s5, s6, s7)
	statistiche = []
	for el1 in tup1:
		stat = []
		for el2 in tup2:
			stat.append(stats.ks_2samp(el1, el2))
		statistiche.append(stat)
	with open('./data/clusterstats/{}.tsv'.format(filename), 'w') as f:
		f.write('\tGroup1({})\tGroup2({})\tGroup3({})\tGroup4({})\tGroup5({})\tGroup6({})\tGroup7({})\n'.format(len(s1), len(s2), len(s3), len(s4), len(s5), len(s6), len(s7)))
		for i, ks in enumerate(statistiche, start = 1):
			f.write('Group{}\t'.format(i))
			line = ''
			for stat in ks:
				line += '{}\t'.format(stat[1])
			line.strip()
			f.write(line + '\n')
		print('Done')


group1 = ward_groups[0]
group2 = ward_groups[1]
group3 = ward_groups[2]
group4 = ward_groups[3]
group5 = ward_groups[4]
group6 = ward_groups[5]
group7 = ward_groups[6]



####################################### Lattato Liquorale #######################################

df = pd.read_excel('../database/db.xlsx', sheet_name = 'SM NM', usecols = 'A,CI')
df.columns = ['patient_id', 'lact']
df1 = df[df['patient_id'].isin(group1)].dropna()
df2 = df[df['patient_id'].isin(group2)].dropna()
df3 = df[df['patient_id'].isin(group3)].dropna()
df4 = df[df['patient_id'].isin(group4)].dropna()
df5 = df[df['patient_id'].isin(group5)].dropna()
df6 = df[df['patient_id'].isin(group6)].dropna()
df7 = df[df['patient_id'].isin(group7)].dropna()
s1 = list(df1['lact'].astype(float))
s2 = list(df2['lact'].astype(float))
s3 = list(df3['lact'].astype(float))
s4 = list(df4['lact'].astype(float))
s5 = list(df5['lact'].astype(float))
s6 = list(df6['lact'].astype(float))
s7 = list(df7['lact'].astype(float))
cross_ks(s1, s2, s3, s4, s5, s6, s7, 'lact')


####################################### Vitamina D #######################################

df = pd.read_excel('../database/db.xlsx', sheet_name = 'SM NM', usecols = 'A,CJ')
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
cross_ks(s1, s2, s3, s4, s5, s6, s7, 'vitd')


####################################### Età #######################################

df = pd.read_excel('../database/db.xlsx', sheet_name = 'SM NM', usecols = 'A,CM')
df.columns = ['patient_id', 'age']
df1 = df[df['patient_id'].isin(group1)].dropna()
df2 = df[df['patient_id'].isin(group2)].dropna()
df3 = df[df['patient_id'].isin(group3)].dropna()
df4 = df[df['patient_id'].isin(group4)].dropna()
df5 = df[df['patient_id'].isin(group5)].dropna()
df6 = df[df['patient_id'].isin(group6)].dropna()
df7 = df[df['patient_id'].isin(group7)].dropna()
s1 = list(df1['age'].astype(float))
s2 = list(df2['age'].astype(float))
s3 = list(df3['age'].astype(float))
s4 = list(df4['age'].astype(float))
s5 = list(df5['age'].astype(float))
s6 = list(df6['age'].astype(float))
s7 = list(df7['age'].astype(float))
cross_ks(s1, s2, s3, s4, s5, s6, s7, 'age')


####################################### EDSS_T0 #######################################

df = pd.read_excel('../database/db.xlsx', sheet_name = 'SM NM', usecols = 'A,DC')
df.columns = ['patient_id', 'edss_t0']
df1 = df[df['patient_id'].isin(group1)].dropna()
df2 = df[df['patient_id'].isin(group2)].dropna()
df3 = df[df['patient_id'].isin(group3)].dropna()
df4 = df[df['patient_id'].isin(group4)].dropna()
df5 = df[df['patient_id'].isin(group5)].dropna()
df6 = df[df['patient_id'].isin(group6)].dropna()
df7 = df[df['patient_id'].isin(group7)].dropna()
s1 = list(df1['edss_t0'].astype(float))
s2 = list(df2['edss_t0'].astype(float))
s3 = list(df3['edss_t0'].astype(float))
s4 = list(df4['edss_t0'].astype(float))
s5 = list(df5['edss_t0'].astype(float))
s6 = list(df6['edss_t0'].astype(float))
s7 = list(df7['edss_t0'].astype(float))
cross_ks(s1, s2, s3, s4, s5, s6, s7, 'edss_t0')


####################################### EDSS_T1 #######################################

df = pd.read_excel('../database/db.xlsx', sheet_name = 'SM NM', usecols = 'A,FF')
df.columns = ['patient_id', 'edss_t1']
df1 = df[df['patient_id'].isin(group1)].dropna()
df2 = df[df['patient_id'].isin(group2)].dropna()
df3 = df[df['patient_id'].isin(group3)].dropna()
df4 = df[df['patient_id'].isin(group4)].dropna()
df5 = df[df['patient_id'].isin(group5)].dropna()
df6 = df[df['patient_id'].isin(group6)].dropna()
df7 = df[df['patient_id'].isin(group7)].dropna()
s1 = list(df1['edss_t1'].astype(float))
s2 = list(df2['edss_t1'].astype(float))
s3 = list(df3['edss_t1'].astype(float))
s4 = list(df4['edss_t1'].astype(float))
s5 = list(df5['edss_t1'].astype(float))
s6 = list(df6['edss_t1'].astype(float))
s7 = list(df7['edss_t1'].astype(float))
cross_ks(s1, s2, s3, s4, s5, s6, s7, 'edss_t1')

####################################### DELTA_EDSS #######################################

df = pd.read_excel('../database/db.xlsx', sheet_name = 'SM NM', usecols = 'A,DC,FF')
df.columns = ['patient_id', 'edss_t0', 'edss_t1']
df1 = df[df['patient_id'].isin(group1)].dropna()
df2 = df[df['patient_id'].isin(group2)].dropna()
df3 = df[df['patient_id'].isin(group3)].dropna()
df4 = df[df['patient_id'].isin(group4)].dropna()
df5 = df[df['patient_id'].isin(group5)].dropna()
df6 = df[df['patient_id'].isin(group6)].dropna()
df7 = df[df['patient_id'].isin(group7)].dropna()
s1a = list(df1['edss_t0'].astype(float))
s2a = list(df2['edss_t0'].astype(float))
s3a = list(df3['edss_t0'].astype(float))
s4a = list(df4['edss_t0'].astype(float))
s5a = list(df5['edss_t0'].astype(float))
s6a = list(df6['edss_t0'].astype(float))
s7a = list(df7['edss_t0'].astype(float))
s1b = list(df1['edss_t1'].astype(float))
s2b = list(df2['edss_t1'].astype(float))
s3b = list(df3['edss_t1'].astype(float))
s4b = list(df4['edss_t1'].astype(float))
s5b = list(df5['edss_t1'].astype(float))
s6b = list(df6['edss_t1'].astype(float))
s7b = list(df7['edss_t1'].astype(float))
def listDiff(l1, l2):
	diff_list = []
	for i in range(len(l1)):
		diff_list.append(l2[i] - l1[i])
	return(diff_list)
s1 = listDiff(s1a, s1b)
s2 = listDiff(s2a, s2b)
s3 = listDiff(s3a, s3b)
s4 = listDiff(s4a, s4b)
s5 = listDiff(s5a, s5b)
s6 = listDiff(s6a, s6b)
s7 = listDiff(s7a, s7b)
cross_ks(s1, s2, s3, s4, s5, s6, s7, 'delta_edss')

####################################### Ricadute alla PL #######################################

df = pd.read_excel('../database/db.xlsx', sheet_name = 'SM NM', usecols = 'A,EZ')
df.columns = ['patient_id', 'rel_t0']
df1 = df[df['patient_id'].isin(group1)].dropna()
df1.replace(to_replace = 'NO', value = 0, inplace = True)
df2 = df[df['patient_id'].isin(group2)].dropna()
df2.replace(to_replace = 'NO', value = 0, inplace = True)
df3 = df[df['patient_id'].isin(group3)].dropna()
df3.replace(to_replace = 'NO', value = 0, inplace = True)
df4 = df[df['patient_id'].isin(group4)].dropna()
df4.replace(to_replace = 'NO', value = 0, inplace = True)
df5 = df[df['patient_id'].isin(group5)].dropna()
df5.replace(to_replace = 'NO', value = 0, inplace = True)
df6 = df[df['patient_id'].isin(group6)].dropna()
df6.replace(to_replace = 'NO', value = 0, inplace = True)
df7 = df[df['patient_id'].isin(group7)].dropna()
df7.replace(to_replace = 'NO', value = 0, inplace = True)
s1 = list(df1['rel_t0'].astype(float))
s2 = list(df2['rel_t0'].astype(float))
s3 = list(df3['rel_t0'].astype(float))
s4 = list(df4['rel_t0'].astype(float))
s5 = list(df5['rel_t0'].astype(float))
s6 = list(df6['rel_t0'].astype(float))
s7 = list(df7['rel_t0'].astype(float))
cross_ks(s1, s2, s3, s4, s5, s6, s7, 'rel_t0')


####################################### MSSS #######################################

df = pd.read_excel('../database/db.xlsx', sheet_name = 'SM NM', usecols = 'A,FG')
df.columns = ['patient_id', 'msss']
df1 = df[df['patient_id'].isin(group1)].dropna()
df2 = df[df['patient_id'].isin(group2)].dropna()
df3 = df[df['patient_id'].isin(group3)].dropna()
df4 = df[df['patient_id'].isin(group4)].dropna()
df5 = df[df['patient_id'].isin(group5)].dropna()
df6 = df[df['patient_id'].isin(group6)].dropna()
df7 = df[df['patient_id'].isin(group7)].dropna()
s1 = list(df1['msss'].astype(float))
s2 = list(df2['msss'].astype(float))
s3 = list(df3['msss'].astype(float))
s4 = list(df4['msss'].astype(float))
s5 = list(df5['msss'].astype(float))
s6 = list(df6['msss'].astype(float))
s7 = list(df7['msss'].astype(float))
cross_ks(s1, s2, s3, s4, s5, s6, s7, 'msss')


####################################### BREMS #######################################

df = pd.read_excel('../database/db.xlsx', sheet_name = 'SM NM', usecols = 'A,JG')
df.columns = ['patient_id', 'brems']
df1 = df[df['patient_id'].isin(group1)].dropna()
df2 = df[df['patient_id'].isin(group2)].dropna()
df3 = df[df['patient_id'].isin(group3)].dropna()
df4 = df[df['patient_id'].isin(group4)].dropna()
df5 = df[df['patient_id'].isin(group5)].dropna()
df6 = df[df['patient_id'].isin(group6)].dropna()
df7 = df[df['patient_id'].isin(group7)].dropna()
s1 = list(df1['brems'].astype(float))
s2 = list(df2['brems'].astype(float))
s3 = list(df3['brems'].astype(float))
s4 = list(df4['brems'].astype(float))
s5 = list(df5['brems'].astype(float))
s6 = list(df6['brems'].astype(float))
s7 = list(df7['brems'].astype(float))
cross_ks(s1, s2, s3, s4, s5, s6, s7, 'brems')


####################################### FSS #######################################

df = pd.read_excel('../database/db.xlsx', sheet_name = 'SM NM', usecols = 'A,QW')
df.columns = ['patient_id', 'fss']
df1 = df[df['patient_id'].isin(group1)].dropna()
df2 = df[df['patient_id'].isin(group2)].dropna()
df3 = df[df['patient_id'].isin(group3)].dropna()
df4 = df[df['patient_id'].isin(group4)].dropna()
df5 = df[df['patient_id'].isin(group5)].dropna()
df6 = df[df['patient_id'].isin(group6)].dropna()
df7 = df[df['patient_id'].isin(group7)].dropna()
s1 = list(df1['fss'].astype(float))
s2 = list(df2['fss'].astype(float))
s3 = list(df3['fss'].astype(float))
s4 = list(df4['fss'].astype(float))
s5 = list(df5['fss'].astype(float))
s6 = list(df6['fss'].astype(float))
s7 = list(df7['fss'].astype(float))
cross_ks(s1, s2, s3, s4, s5, s6, s7, 'fss')


####################################### BMI #######################################

df = pd.read_excel('../database/db.xlsx', sheet_name = 'SM NM', usecols = 'A,UR')
df.columns = ['patient_id', 'bmi']
df1 = df[df['patient_id'].isin(group1)].dropna()
df2 = df[df['patient_id'].isin(group2)].dropna()
df3 = df[df['patient_id'].isin(group3)].dropna()
df4 = df[df['patient_id'].isin(group4)].dropna()
df5 = df[df['patient_id'].isin(group5)].dropna()
df6 = df[df['patient_id'].isin(group6)].dropna()
df7 = df[df['patient_id'].isin(group7)].dropna()
s1 = list(df1['bmi'].astype(float))
s2 = list(df2['bmi'].astype(float))
s3 = list(df3['bmi'].astype(float))
s4 = list(df4['bmi'].astype(float))
s5 = list(df5['bmi'].astype(float))
s6 = list(df6['bmi'].astype(float))
s7 = list(df7['bmi'].astype(float))
cross_ks(s1, s2, s3, s4, s5, s6, s7, 'bmi')


####################################### Ultima ricaduta #######################################

df = pd.read_excel('../database/db.xlsx', sheet_name = 'SM NM', usecols = 'A,CS')
df.columns = ['patient_id', 'last_rel']
df1 = df[df['patient_id'].isin(group1)].dropna()
df2 = df[df['patient_id'].isin(group2)].dropna()
df3 = df[df['patient_id'].isin(group3)].dropna()
df4 = df[df['patient_id'].isin(group4)].dropna()
df5 = df[df['patient_id'].isin(group5)].dropna()
df6 = df[df['patient_id'].isin(group6)].dropna()
df7 = df[df['patient_id'].isin(group7)].dropna()
s1 = list(df1['last_rel'].astype(float))
s2 = list(df2['last_rel'].astype(float))
s3 = list(df3['last_rel'].astype(float))
s4 = list(df4['last_rel'].astype(float))
s5 = list(df5['last_rel'].astype(float))
s6 = list(df6['last_rel'].astype(float))
s7 = list(df7['last_rel'].astype(float))
cross_ks(s1, s2, s3, s4, s5, s6, s7, 'last_rel')

####################################### Volume Lesioni #######################################

df = pd.read_excel('../database/db.xlsx', sheet_name = 'SM NM', usecols = 'A,FQ')
df.columns = ['patient_id', 'les_vol']
df1 = df[df['patient_id'].isin(group1)].dropna()
df2 = df[df['patient_id'].isin(group2)].dropna()
df3 = df[df['patient_id'].isin(group3)].dropna()
df4 = df[df['patient_id'].isin(group4)].dropna()
df5 = df[df['patient_id'].isin(group5)].dropna()
df6 = df[df['patient_id'].isin(group6)].dropna()
df7 = df[df['patient_id'].isin(group7)].dropna()
s1 = list(df1['les_vol'].astype(float))
s2 = list(df2['les_vol'].astype(float))
s3 = list(df3['les_vol'].astype(float))
s4 = list(df4['les_vol'].astype(float))
s5 = list(df5['les_vol'].astype(float))
s6 = list(df6['les_vol'].astype(float))
s7 = list(df7['les_vol'].astype(float))
cross_ks(s1, s2, s3, s4, s5, s6, s7, 'les_vol')