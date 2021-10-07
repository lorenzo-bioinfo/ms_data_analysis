#testing differences in snps alleles frequencies between
#controls and european population using pearson's
#chi-squared test.

import pandas as pd
from scipy import stats

#reading list of available snps

with open('./data/snp/snp_list.csv', 'r') as f:
	snp_list = f.readline().strip().split(',')

#reading table containing snps genotype frequencies
#in european population

snp_df = pd.read_csv('./data/snp/snp_hw.tsv', sep = '\t', header = 0, index_col = 0)

#for each snp, extracting frequencies of the two possible alleles

euro_gen_freq = {} #frequencies I will use as expected in tests

for snp in snp_list:
	p = list(snp_df[snp_df['id'] == snp]['p'])[0]
	p_freq = list(snp_df[snp_df['id'] == snp]['freq_p'])[0]
	q = list(snp_df[snp_df['id'] == snp]['q'])[0]
	q_freq = list(snp_df[snp_df['id'] == snp]['freq_q'])[0]
	euro_gen_freq[snp] = [(p, p_freq), (q, q_freq)]

for snp in euro_gen_freq:
	print(f'{snp} --> {euro_gen_freq[snp]}')

#getting patients data

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

ctrl_df = df[(df['class_int'] == 6) | (df['class_int'] == 7)]
print(ctrl_df)

#now every snp must to be turned into a string,
#made by the concatenation of the genotypes of each
#single patient

snp_strings = {}
for snp in snp_list:
	snp_col = list(ctrl_df[snp].dropna())
	snp_str = ''.join(snp_col)
	snp_strings[snp] = snp_str

#now generating data for chi_squared test

snp_list = list(snp_strings.keys())
chisq_data = {}
for snp in snp_list:
	p = euro_gen_freq[snp][0][0]
	q = euro_gen_freq[snp][1][0]
	p_eurofreq = euro_gen_freq[snp][0][1]
	q_eurofreq = euro_gen_freq[snp][1][1]
	p_obs = snp_strings[snp].count(p)
	q_obs = snp_strings[snp].count(q)
	p_exp = p_eurofreq * len(snp_strings[snp])
	q_exp = q_eurofreq * len(snp_strings[snp])
	if(p_obs != 0) or (q_obs != 0):
		chisq_data[snp] = ([p_obs, q_obs], [p_exp, q_exp], len(snp_strings[snp]))
snps_chisq = list(chisq_data.keys())
print('################# Data check ###############')
for key in chisq_data:
	print(f'{key} --> {chisq_data[key]}')
#now testing

f = open('./data/snp/tests/alleles_ctrl_vs_euro.tsv', 'w')
f.write('snp_id\tsamp_dim\tChi-stat\tp-val\n')

for snp in snps_chisq:
	stat = stats.chisquare(f_obs = chisq_data[snp][0], f_exp = chisq_data[snp][1])
	f.write('{}\t{}\t{}\t{}\n'.format(snp, chisq_data[snp][2], stat[0], stat[1]))
f.close()

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

df = pd.read_csv('./data/snp/tests/alleles_ctrl_vs_euro.tsv', sep = '\t')
df_corr = benjamini_hochberg(df, 'p-val')
df_corr.to_csv('./data/snp/tests/alleles_ctrl_vs_euro_bh.tsv', sep = '\t')
print(df_corr)