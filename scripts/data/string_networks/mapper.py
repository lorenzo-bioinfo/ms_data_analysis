import requests
from time import sleep

url_template = 'https://string-db.org/api/tsv/get_string_ids?identifiers=[{}]&species=9606&echo_query=1&limit=1'

with open('uniprot.txt', 'r') as f:
	uni_list = []
	for line in f:
		uni_list.append(line.strip().split(',')[1])


for uni in uni_list:
	url = url_template.format(uni)
	print('Getting mapping for {}'.format(uni))
	handle = requests.get(url).text
	print('Done:')
	print(handle)
	print('Waiting two secs...')
	sleep(2)