import pandas as pd
import re
import requests
from time import sleep
#importing snps(void columns are omitted)
df = pd.read_excel('../database/db.xlsx', sheet_name = 'SM NM', usecols = 'A,JI:PW')
print(list(df.columns)[0])
#dropping control who is messing up all over the place :¬|
df = df[df['Codice PLASMA NM'] != 'C16']
#now dropping patient ids column
df = df.drop('Codice PLASMA NM', axis = 1)
#dropping empty rows
df.dropna(axis = 0, how = 'all', inplace = True)
#dropping empty columns
df.dropna(axis = 1, how = 'all', inplace = True)
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
log = open('./data/snp/download_error_log.txt', 'w')
f.write('id\tp\tfreq_p\tq\tfreq_q\n')
for snp in snp_list:
	#formatting url to retrieve information about snp
	url = 'https://www.ncbi.nlm.nih.gov/snp/{}'.format(snp)
	#downloading html page with snp informations
	print('Downloading {} page'.format(snp))
	page = requests.get(url).text
	print('OK')
	#getting only table with population frequencies
	try:
		table = pd.read_html(page, match = 'Europe')[0]
		table.drop(0, axis = 0, inplace = True) #doesn't work without this
		#selecting row with allele/frequency in european population
		serie = table[(table['Population'] == 'European') | (table['Population'] == 'Europe')]
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
	except ValueError:
		print(f'Errors with {snp}, adding url to error log')
		log.write(f'{url}\n')
		p = '?'
		freq_p = '?'
		q = '?'
		freq_q = '?'
	f.write('{}\t{}\t{}\t{}\t{}\n'.format(snp, p, freq_p, q, freq_q))
	print('Finished with {}, taking a break...'.format(snp))
	sleep(1)
f.close()
log.close()

'''
df = pd.read_csv('./data/snp/snp_info.tsv', sep = '\t')
print(df)