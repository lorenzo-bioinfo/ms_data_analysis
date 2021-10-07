import pandas as pd
from time import sleep
import requests

def get_string_network(str_id, score, nodes, dictionary):
	url_model = 'https://string-db.org/api/tsv/network?identifiers={}&species=9606&required_score={}&add_nodes={}'
	url = url_model.format(str_id, score, nodes)
	handle = requests.get(url)
	with open('str_{}.tsv'.format(dictionary[str_id]), 'w') as f:
		f.write(handle.text)

with open('str_list.txt', 'r') as f:
	cyt_list = f.readline().split(',')

with open('uniprot.txt', 'r') as f:
	str_dict = {}
	str_list = []
	for line in f:
		cyt = line.strip().split(',')[0]
		string = line.strip().split(',')[2].split('.')[1]
		str_list.append(string)
		str_dict[string] = cyt

for strid in str_list:
	print('Downloading {} network'.format(str_dict[strid]))
	get_string_network(strid, 200, 500, str_dict)
	print('Done, taking a break..')
	sleep(2)