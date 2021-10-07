import pandas as pd
import random
from scipy.cluster import hierarchy
from matplotlib import pyplot as plt
import seaborn as sns

import pandas as pd
from scipy import stats

#reading list of available snps

with open('./data/snp/snp_list.csv', 'r') as f:
	snp_list = f.readline().strip().split(',')

#reading table containing snps genotype frequencies
#in european population

snp_df = pd.read_csv('./data/snp/snp_hw.tsv', sep = '\t', header = 0, index_col = 0)

#for each snp, extracting frequencies of 3 possible genotypes

euro_gen_freq = {} #frequencies I will use as expected in tests

for snp in snp_list:
	pp = list(snp_df[snp_df['id'] == snp]['pp'])[0]
	pq = list(snp_df[snp_df['id'] == snp]['pq'])[0]
	qq = list(snp_df[snp_df['id'] == snp]['qq'])[0]
	euro_gen_freq[snp] = [pp, pq, qq]

#WHERE:
#	pp is equal to omozygote for most frequent allele
#	pq is eterozygote
#	qq is omozygote for least frequent allele

#now I will create a list of all 3 possible genotypes for each snp
	#this will contains the 2 alleles for each snp
snp_alleles = []
for snp in snp_list:
	entry = []
	entry.append(snp)
	entry.append(list(snp_df[snp_df['id'] == snp]['p'])[0])
	entry.append(list(snp_df[snp_df['id'] == snp]['q'])[0])
	snp_alleles.append(tuple(entry))
	#this dictionary will contain the possible genotypes.
	#since I cannot know wether the eterozygotes are encoded
	#as pq or qp, I'll have to add both the formattings
snp_genotp = {}
for entry in snp_alleles:
	omop = entry[1] + entry[1]
	eth1 = entry[1] + entry[2]
	eth2 = entry[2] + entry[1]
	omoq = entry[2] + entry[2]
	snp_genotp[entry[0]] = (omop, eth1, eth2, omoq)


#now I'm going to read the data about snps from the 'database'

df = pd.read_excel('../database/db.xlsx', sheet_name = 'SM NM', usecols = 'A,CB,JI:PW')
#dropping control who is messing up all over the place :¬|
df = df[df['Codice PLASMA NM'] != 'C16']
#dropping empty columns
df.dropna(axis = 1, how = 'all', inplace = True)
#renaming columns because comfort
snp_list.insert(0, 'patient_id')
snp_list.insert(1, 'class_int')
df.columns = snp_list
snp_list.pop(0)
snp_list.pop(0)

#Isolating RR patients
rr_df = df[df['class_int'] == 3]
#rr_df.index = list(rr_df['patient_id'])
#rr_df.drop('patient_id', axis = 1, inplace = True)
#rr_df = rr_df.dropna(axis = 0, how = 'all')
#rr_df.fillna('Empty', inplace = True)

#selecting only snps which showed significant differences
#between rr and european population, but not between eu and ctrls

snps_rr = 'patient_id,rs11913721,rs1474347,rs11568817,rs1554973,rs1799880,rs2069763,rs3757351,rs2227306,rs2069772,rs7301328,rs7598440,rs743409,rs4680,rs12722489,rs2227284,rs1800797,rs2283792,rs4251961,rs10767664,rs2069812,rs1927914,rs7131056'.split(',')
df_ok = rr_df[rr_df.columns.intersection(snps_rr)].dropna()
snps_rr.remove('patient_id')
df_ok.index = df_ok['patient_id']
df_ok.drop('patient_id', axis = 1, inplace = True)
df_ok = df_ok.dropna(axis = 0, how = 'all')
df_ok['patient_id'] = df_ok.index
pat_list = list(df_ok['patient_id'])
print(df_ok, pat_list)

#so there are some nans  and I want to fill them.
#The best idea I've had so far is to fill nans with random
#genotypes, based on frequencies observed in patients who have
#a value for that snp

rr_freqs = {}
for snp in snps_rr:
	arr = list(df_ok[snp])
	arr_clean = [x for x in arr if x != 'Empty']
	genotps = set(arr_clean)
	snps = {}
	for el in genotps:
		snps[el] = arr.count(el) / len(arr_clean)
	rr_freqs[snp] = snps

#creating samples from which I'll extract random genotypes
#for each snp

rr_samples = {}
for snp in rr_freqs:
	sample = []
	for gentp in rr_freqs[snp]:
		dim = int(rr_freqs[snp][gentp] * 100)
		lista = (f'{gentp},' * dim).split(',')
		lista_ok = [x for x in lista if x != '']
		sample.extend(lista_ok)
	rr_samples[snp] = sample

#recoding genotypes for clustering:
	#pp genotypes will be 0,
	#pq genotypes will be 1,
	#qq genotypes will be 2

recoded = {}
for pat in pat_list:
	row = df_ok[df_ok['patient_id'] == pat]
	pat_rec = []
	for snp in snps_rr:
		omop = snp_genotp[snp][0]
		eth1 = snp_genotp[snp][1]
		eth2 = snp_genotp[snp][2]
		omoq = snp_genotp[snp][3]
		if row[snp].values[0] == omop:
			pat_rec.append(0)
		elif (row[snp].values[0] == eth1) or (row[snp].values[0] == eth2):
			pat_rec.append(1)
		elif row[snp].values[0] == omoq:
			pat_rec.append(2)
		else:
			random_gen = random.choice(rr_samples[snp])
			random_rev = random_gen[::-1]
			print(random_gen)
			if random_gen == omop:
				pat_rec.append(0)
			elif random_gen == omoq:
				pat_rec.append(2)
			else:
				pat_rec.append(1)
	recoded[pat] = pat_rec

df = pd.DataFrame.from_dict(recoded, orient = 'index', columns = snps_rr)

#now clustering

cluster_col = hierarchy.linkage(df.T, method="ward", metric="euclidean")
cluster_row = hierarchy.linkage(df, method="ward", metric="euclidean")
clusterfig = sns.clustermap(df, row_linkage = cluster_row, col_linkage = cluster_col, cmap = 'mako')
index_col = clusterfig.dendrogram_col.reordered_ind
index_row = clusterfig.dendrogram_row.reordered_ind
plt.savefig('../plots/snp_clusters/rr_only_vs_euro_per_slide.png', dpi = 300)
plt.clf()
clusterfig_full = sns.clustermap(df, row_linkage = cluster_row, col_linkage = cluster_col, cmap = 'mako', yticklabels = True, figsize = (10, len(df)/4))
plt.savefig('../plots/snp_clusters/rr_only_vs_euro.png', dpi = 300)
plt.clf()

#getting list of patient ids (not ordered)
ids_list = list(df.index)
ids_ordered = []
#using the indexes provided by index_row to obtain
#a list of ordered ids (as in the cluster figure)
print('INDEX ROW CHECK')
print(index_row)
for i in index_row:
	ids_ordered.append(ids_list[i])

#exporting ordered patients ids

with open('data/cluster_groups/ordered_rr_snps.txt', 'w') as f:
	for idn in ids_ordered:
		f.write(str(idn) + '\n')