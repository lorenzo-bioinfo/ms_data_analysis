#testing differences in snps genotypes between
#RR patients and european population using pearson's
#chi-squared test.

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

#adding a check, just to see if this is running :)
for entry in snp_genotp:
	print(f'{entry} --> {snp_genotp[entry]}')

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
#checking again
print(df)

#Isolating controls

rr_df = df[df['class_int'] == 3]
print(rr_df)

#now I will proceed to recode the genotypes using integers.
#In particular, pp genotype will be 0, pq 1 and qq 2.
#This will be useful for testing, clustering and other stuff.

df_samp = {}
for entry in snp_genotp:
	samp = []
	snp_row = list(rr_df[entry].dropna())
	recode = []
	for gen in snp_row:
		if gen == snp_genotp[entry][0]:
			recode.append(0)
		elif gen == snp_genotp[entry][1] or gen == snp_genotp[entry][2]:
			recode.append(1)
		elif gen == snp_genotp[entry][3]:
			recode.append(2)
	df_samp[entry] = recode

#turns out the genotypes for some of the snps in our db
#do not match the 'possible' alleles reported on ncbi.
#Since I trust NCBI way more than our data, I will ignore
#those snps and move on only with the 'good' ones

clean_samp = dict(df_samp)

for snp in df_samp:
	if len(df_samp[snp]) == 0:
		del clean_samp[snp]
for snp in clean_samp:
	print(snp, clean_samp[snp], len(clean_samp[snp]))

#and now getting a list of the remaining snps

snps_clean = list(clean_samp.keys())
#saving clean snps list to file
with open('./data/snp/snp_list_clean.csv', 'w') as f:
	f.write(','.join(snps_clean))

#now I will get, for each snp, the observed frequencies (by counting
#how many times each genotype is present in controls) and the expected
#frequencies, by multiplyin the genotype expected relative frequency
#in the european population by the length of the controls' sample

samp_gen_freq = {}

for snp in snps_clean:
	tot = len(clean_samp[snp])
	pp = clean_samp[snp].count(0)
	pq = clean_samp[snp].count(1)
	qq = clean_samp[snp].count(2)
	pp_exp = euro_gen_freq[snp][0]*tot
	pq_exp = euro_gen_freq[snp][1]*tot
	qq_exp = euro_gen_freq[snp][2]*tot
	samp_gen_freq[snp] = [[pp, pq, qq], [pp_exp, pq_exp, qq_exp], tot]
	print(samp_gen_freq[snp])

#and now finally testing

f = open('./data/snp/tests/rr_vs_euro.tsv', 'w')
f.write('snp_id\tsamp_dim\tChi-stat\tp-val\n')

for snp in snps_clean:
	stat = stats.chisquare(f_obs = samp_gen_freq[snp][0], f_exp = samp_gen_freq[snp][1])
	f.write('{}\t{}\t{}\t{}\n'.format(snp, samp_gen_freq[snp][2], stat[0], stat[1]))
f.close()

#applying benjamini-hochberg correction to p-values

def benjamini_hochberg(dataframe, dataframe_column):
	#dataframe_column is the name or index of the col containing p-values
	dataframe.sort_values(by = dataframe_column, inplace = True)
	p_values = list(dataframe[dataframe_column].astype(float))
	q_values = []
	n = len(p_values)
	for i, p_value in enumerate(p_values):
		q_values.append((p_value * n)/(i+1))
	#adding q_values column to dataframe
	dataframe['q_values'] = q_values
	dataframe.sort_values('q_values', inplace = True)
	return dataframe

df = pd.read_csv('./data/snp/tests/rr_vs_euro.tsv', sep = '\t')
df_corr = benjamini_hochberg(df, 'p-val')
df_corr.to_csv('./data/snp/tests/rr_vs_euro_bh.tsv', sep = '\t')
print(df_corr)