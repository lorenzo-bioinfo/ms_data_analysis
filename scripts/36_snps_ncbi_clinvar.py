import requests
import pandas as pd
from bs4 import BeautifulSoup as bs

#getting list of snps ids from file
with open('./data/snp/snp_list.csv', 'r') as f:
	snp_list = f.readline().strip().split(',')

f = open('./data/snp/snp_clinvar.tsv', 'w')
f.write('snp_id\treported\tclin_significance\tgene_1\tgene_consequence1\tgene_2\tgene_consequence2\n')
for snp in snp_list:
	#formatting url to retrieve information about snp
	url = 'https://www.ncbi.nlm.nih.gov/snp/{}'.format(snp)
	#downloading html page with snp informations
	print('Downloading {} page'.format(snp))
	page = requests.get(url).text
	print('OK')
	#using bs4 to parse the document
	soup = bs(page, 'html.parser') #take the whole page
	#finding section containing clinical information:
	try:
		area = soup.find_all('dl', class_='usa-width-one-half')[1]
		#extracting significance and consequence:
		try:
			reported = area.find('dd').text.strip()
		except AttributeError:
			reported = None
		try:
			info = area.find('span').text.strip()
			gene1 = info.split(':')[0]
			consequence1 = info.split(':')[1]
		except AttributeError:
			consequence1 = None
			gene1 = None
		except IndexError:
			consequence1 = None
			gene1 = None
		try:
			info = area.find('div').text.strip()
			consequence2 = info.split(':')[1]
			gene2 = info.split(':')[0]
		except AttributeError:
			consequence2 = None
			gene2 = None
		except IndexError:
			consequence1 = None
			gene1 = None
		if reported == 'Reported in ClinVar':
			table = soup.find('table', id = "clinical_significance_datatable")
			lista = table.find_all('td')
			significance = str(lista[len(lista) - 1])
		else:
			significance = 'Not reported'
	except IndexError:
		reported = None
		significance = 'Not reported'
		gene1 = None
		consequence1 = None
		gene2 = None
		consequence2 = None

	f.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(snp, reported, significance, gene1, consequence1, gene2, consequence2))
f.close()