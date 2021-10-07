#summarizing the results of snps testing

import pandas as pd

all_ctrl_euro = pd.read_csv('./data/snp/tests/alleles_ctrl_vs_euro_bh.tsv', sep = '\t', index_col = 0, header = 0)
all_rr_euro = pd.read_csv('./data/snp/tests/alleles_rr_vs_euro_bh.tsv', sep = '\t', index_col = 0, header = 0)
all_rr_ctrl = pd.read_csv('./data/snp/tests/alleles_rr_vs_ctrl_bh.tsv', sep = '\t', index_col = 0, header = 0)
ctrl_euro = pd.read_csv('./data/snp/tests/ctrl_vs_euro_bh.tsv', sep = '\t', index_col = 0, header = 0)
rr_euro = pd.read_csv('./data/snp/tests/rr_vs_euro_bh.tsv', sep = '\t', index_col = 0, header = 0)
rr_ctrl = pd.read_csv('./data/snp/tests/ctrl_vs_rr_bh.tsv', sep = '\t', index_col = 0, header = 0)

all_ctrl_euro_sign = set(list(all_ctrl_euro[all_ctrl_euro['q_values'] < 0.05]['snp_id']))
all_rr_euro_sign = set(list(all_rr_euro[all_rr_euro['q_values'] < 0.05]['snp_id']))
all_rr_ctrl_sign = set(list(all_rr_ctrl[all_rr_ctrl['q_values'] < 0.05]['snp_id']))
ctrl_euro_sign = set(list(ctrl_euro[ctrl_euro['q_values'] < 0.05]['snp_id']))
rr_euro_sign = set(list(rr_euro[rr_euro['q_values'] < 0.05]['snp_id']))
rr_ctrl_sign = set(list(rr_ctrl[rr_ctrl['q_values'] < 0.05]['snp_id']))

list_a = [all_ctrl_euro_sign, all_rr_euro_sign, all_rr_ctrl_sign, ctrl_euro_sign, rr_euro_sign, rr_ctrl_sign]
list_b = [all_ctrl_euro_sign, all_rr_euro_sign, all_rr_ctrl_sign, ctrl_euro_sign, rr_euro_sign, rr_ctrl_sign]
matrix = []
for i in range(len(list_a)):
	row = []
	for j in range(len(list_b)):
		row.append(len(list_a[i] & list_b[j]))
	matrix.append(row)
for row in matrix:
	print(row)

indexes = ['a_ct_eu ({})', 'a_rr_eu ({})', 'a_ct_rr ({})', 'ct_eu ({})', 'rr_eu ({})', 'ct_rr ({})']
indexes_ok = []
for i in range(len(indexes)):
	indexes_ok.append(indexes[i].format(len(list_a[i])))

mat = pd.DataFrame(matrix)
mat.index = indexes_ok
mat.columns = indexes_ok
print(mat)
mat.to_csv('./data/snp/tests/summary/tests_summary_table.tsv', sep = '\t')

with open('./data/snp/tests/summary/tests_summary_snps.txt', 'w') as f:
	for i in range(len(indexes_ok)):
		f.write(f'# {indexes_ok[i]}\n')
		for snp in list_a[i]:
			f.write(f'{snp}\n')

with open('./data/snp/tests/summary/tests_summary_common.txt', 'w') as f:
	for i in range(len(indexes_ok)):
		for j in range(len(indexes_ok)):
			com = list(list_a[i] & list_b[j])
			f.write(f'# {indexes_ok[i]} |{len(com)}| {indexes_ok[j]}\n')
			f.write('\n'.join(com) + '\n')

gwas = pd.read_csv('./data/snp/gwas_info.tsv', sep = '\t')

snp_list = list(gwas['SNPS'].astype(str))
gwas_diseases = {}
for snp in snp_list:
	disease_list = list(gwas[gwas['SNPS'] == str(snp)]['DISEASE/TRAIT'])
	gwas_diseases[snp] = disease_list

with open('./data/snp/tests/summary/tests_summary_gwas.txt', 'w') as f:
	for i in range(len(indexes_ok)):
		f.write(f'# {indexes_ok[i]}\n')
		for snp in list_a[i]:
			try:
				f.write(f'{snp} - {"|".join(gwas_diseases[snp])}\n')
			except KeyError:
				f.write(f'{snp} - Not in GWAS\n')

with open('./data/snp/tests/summary/eu_rr_only.txt', 'w') as f:
	rr_only = [snp for snp in all_rr_euro_sign if snp not in all_ctrl_euro_sign]
	f.write(f'#Alleles({len(rr_only)})\n')
	for snp in rr_only:
		f.write(f'{snp}\n')
	rr_genotp = [snp for snp in rr_euro_sign if snp not in ctrl_euro_sign]
	f.write(f'#Genotypes({len(rr_genotp)})\n')
	for snp in rr_genotp:
		f.write(f'{snp}\n')
	rr_common = set(rr_only) & set(rr_genotp)
	f.write(f'#Common({len(rr_common)})\n')
	for snp in rr_common:
		f.write(f'{snp}\n')

with open('./data/snp/tests/summary/rr_only_gwas.txt', 'w') as f:
	for snp in rr_common:
		try:
			f.write(f'{snp} - {"|".join(gwas_diseases[snp])}\n')
		except KeyError:
			f.write(f'{snp} - Not in GWAS\n')

clinvar_df = pd.read_csv('./data/snp/snp_clinvar.tsv', sep = '\t')
common_df = clinvar_df[clinvar_df['snp_id'].isin(rr_common)]
print(common_df)
common_df.to_csv('./data/snp/tests/summary/rr_only_clinvar.tsv', sep = '\t')