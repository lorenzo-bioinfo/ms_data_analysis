import requests
import pandas as pd
import re
#importing snps(void columns are omitted)
df = pd.read_excel('../database/db.xlsx', sheet_name = 'SM NM', usecols = 'LN:LV, LY, LZ, MC, MH, MN, MQ, MT, MZ, NA, NC, NF, NG, NL, NM, NO, NP, NV, NY, OC, OD, OI, OJ, OM, ON, OS:OU, OW:SB')
#dropping void rows
df.dropna(axis = 0, how = 'all', inplace = True)
#getting list of snps ids and saving them into a file for fast recovery
dirty_list = list(df.columns)
snp_list = []
for entry in dirty_list:
	snp_list.append((re.findall(r'rs\d+', entry)[0]))
print(snp_list)
with open('./data/snp/snp_list.csv', 'w') as f:
	f.write(','.join(snp_list))

###################################################################################
# WARNING: the following code has problems due to the format of alternative		  #
# alleles on NCBI page. The file generated has been corrected by hand. DON'T run  #
# the script again before making a copy of snp_info.tsv                           #
###################################################################################

'''
f = open('./data/snp/snp_info.tsv', 'w')
f.write('id\tp\tfreq_p\tq\tfreq_q\n')
for snp in snp_list:
	#formatting url to retrieve information about snp
	url = 'https://www.ncbi.nlm.nih.gov/snp/{}'.format(snp)
	#downloading html page with snp informations
	print('Downloading {} page'.format(snp))
	page = requests.get(url).text
	print('OK')
	#getting only table with population frequencies
	table = pd.read_html(page, match = 'European')[0]
	table.drop(0, axis = 0, inplace = True) #doesn't work without this
	#selecting row with allele/frequency in european population
	serie = table[table['Population'] == 'European']
	#extracting allele and frequency information
	ref, ref_freq = list(serie['Ref Allele'])[0].split('=')[0], list(serie['Ref Allele'])[0].split('=')[1]
	alt, alt_freq = list(serie['Alt Allele'])[0].split('=')[0], list(serie['Alt Allele'])[0].split('=')[1]
	if ref_freq > alt_freq:
		p = ref
		q = alt
		freq_p = ref_freq
		freq_q = alt_freq
	else:
		p = alt
		q = ref
		freq_p = alt_freq
		freq_q = ref_freq
	f.write('{}\t{}\t{}\t{}\t{}\n'.format(snp, p, freq_p, q, freq_q))
f.close()
'''

df = pd.read_csv('./data/snp/snp_info.tsv', sep = '\t')
print(df)