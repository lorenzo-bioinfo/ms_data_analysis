import pandas as pd

all_ctrl_euro = pd.read_csv('./../data/snp/tests/alleles_ctrl_vs_euro_bh.tsv', sep = '\t', index_col = 0, header = 0)
all_rr_euro = pd.read_csv('./../data/snp/tests/alleles_rr_vs_euro_bh.tsv', sep = '\t', index_col = 0, header = 0)
all_rr_ctrl = pd.read_csv('./../data/snp/tests/alleles_rr_vs_ctrl_bh.tsv', sep = '\t', index_col = 0, header = 0)
ctrl_euro = pd.read_csv('./../data/snp/tests/ctrl_vs_euro_bh.tsv', sep = '\t', index_col = 0, header = 0)
rr_euro = pd.read_csv('./../data/snp/tests/rr_vs_euro_bh.tsv', sep = '\t', index_col = 0, header = 0)
rr_ctrl = pd.read_csv('./../data/snp/tests/ctrl_vs_rr_bh.tsv', sep = '\t', index_col = 0, header = 0)

all_ctrl_euro_sign = set(list(all_ctrl_euro[all_ctrl_euro['q_values'] < 0.05]['snp_id']))
all_rr_euro_sign = set(list(all_rr_euro[all_rr_euro['q_values'] < 0.05]['snp_id']))
all_rr_ctrl_sign = set(list(all_rr_ctrl[all_rr_ctrl['q_values'] < 0.05]['snp_id']))
ctrl_euro_sign = set(list(ctrl_euro[ctrl_euro['q_values'] < 0.05]['snp_id']))
rr_euro_sign = set(list(rr_euro[rr_euro['q_values'] < 0.05]['snp_id']))
rr_ctrl_sign = set(list(rr_ctrl[rr_ctrl['q_values'] < 0.05]['snp_id']))

#all_rr_euro_sign
#all_ctrl_euro_sign
#rr_euro_sign
#ctrl_euro_sign
#rr_ctrl_sign

rr_euro_ok = all_rr_euro_sign & rr_euro_sign
rr_ct_ok = rr_ctrl_sign & all_ctrl_euro_sign
ct_euro_ok = all_ctrl_euro_sign & ctrl_euro_sign

gwas = pd.read_csv('./../data/snp/gwas_info.tsv', sep = '\t')

snp_list = list(gwas['SNPS'].astype(str))
gwas_diseases = {}
for snp in snp_list:
	disease_list = list(gwas[gwas['SNPS'] == str(snp)]['DISEASE/TRAIT'])
	gwas_diseases[snp] = disease_list

#gwas disease for rr vs euro sign:
f = open('snps_diseases.txt', 'w')
f.write('############# RR v EU #############\n')
rr_euro_disease = []
for snp in rr_euro_ok:
	try:
		rr_euro_disease.extend(gwas_diseases[snp])
	except KeyError:
		pass
rr_euro_disease = set(rr_euro_disease)		
for d in rr_euro_disease:
	f.write(f'{d}\n')

f.write('############# RR v CT #############\n')
rr_ct_disease = []
for snp in rr_ct_ok:
	try:
		rr_ct_disease.extend(gwas_diseases[snp])
	except KeyError:
		pass
rr_ct_disease = set(rr_ct_disease)
for d in rr_ct_disease:
	f.write(f'{d}\n')

f.write('############# CT v EU #############\n')
ct_eu_disease = []
for snp in set(ct_euro_ok):
	try:
		ct_eu_disease.extend(gwas_diseases[snp])
	except KeyError:
		pass
ct_eu_disease = set(ct_eu_disease)
for d in ct_eu_disease:
	f.write(f'{d}\n')
f.close()