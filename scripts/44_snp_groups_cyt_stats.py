import pandas as pd
from scipy import stats
from numpy import median

#defining function to get groups

def getGroups(file, stops):
	ordered_ids = []
	with open(file) as f:
		for line in f:
			ordered_ids.append(int(line.strip()))
	mask = []
	for i in range(len(ordered_ids)):
		if ordered_ids[i] in stops:
			mask.append(-1)
		else:
			mask.append(0)
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


#extracting groups from snps genotypes based clustering (ward linkage method)
snp_stops = [862, 826]
snp_groups = getGroups('./data/cluster_groups/ordered_rr_snps.txt', snp_stops)
for gp in snp_groups:
	print(len(gp))

def cross_ks(s1, s2, s3, filename):
	tup1 = (s1, s2, s3)
	tup2 = (s1, s2, s3)
	statistiche = []
	for el1 in tup1:
		stat = []
		for el2 in tup2:
			stat.append(stats.ks_2samp(el1, el2))
		statistiche.append(stat)
	with open('./data/clusterstats_snps_cyt/{}.tsv'.format(filename), 'w') as f:
		f.write('\tGroup1({})\tGroup2({})\tGroup3({})\n'.format(len(s1), len(s2), len(s3)))
		for i, ks in enumerate(statistiche, start = 1):
			f.write('Group{}\t'.format(i))
			line = ''
			for stat in ks:
				line += '{}\t'.format(stat[1])
			line.strip()
			f.write(line + '\n')


group1 = snp_groups[0]
group2 = snp_groups[1]
group3 = snp_groups[2]

cyt_list = 'IL1B,IL2,IL4,IL5,IL6,IL7,CXCL8,IL10,IL12B,IL13,IL17A,CSF3,CSF2,IFNG,CCL2,CCL4,TNF,IL1RN,IL9,IL15,CCL11,FGF2,CXCL10,PDGFB,CCL5,VEGFA,CCL3'.split(',')
cyt_cols = 'F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,AA,AB,AC,AD,AE,AF'.split(',')

for i in range(len(cyt_list)):
	cols_use = f'A,{cyt_cols[i]}'
	cyt = cyt_list[i]
	df = pd.read_excel('../database/db.xlsx', sheet_name = 'SM NM', usecols = cols_use)
	df.columns = ['patient_id', cyt]
	df1 = df[df['patient_id'].isin(group1)].dropna()
	df2 = df[df['patient_id'].isin(group2)].dropna()
	df3 = df[df['patient_id'].isin(group3)].dropna()
	s1 = list(df1[cyt].astype(float))
	s2 = list(df2[cyt].astype(float))
	s3 = list(df3[cyt].astype(float))
	cross_ks(s1, s2, s3, cyt)