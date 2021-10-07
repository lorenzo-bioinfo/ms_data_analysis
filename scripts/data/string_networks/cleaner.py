### add case for cytokin in second node

import pandas as pd

prova_df = pd.read_csv('str_TNF.tsv', sep = '\t')
print(prova_df.head(20))

def extract_spec_interac(cyt):
	full_int_df = pd.read_csv('str_{}.tsv'.format(cyt), sep = '\t')
	min_score = full_int_df['score'].min()
	int_dict = {}
	df1_dirt = full_int_df[full_int_df['preferredName_A'] == cyt]
	df2_dirt = full_int_df[full_int_df['preferredName_B'] == cyt]
	minim1 = df1_dirt['score'].min()
	minim2 = df2_dirt['score'].min()
	df1 = df1_dirt[df1_dirt['score'] > minim1]
	df2 = df2_dirt[df2_dirt['score'] > minim2]
	list1 = list(df1['preferredName_B'])
	list2 = list(df2['preferredName_A'])
	score1 = list(df1['score'].astype(float))
	score2 = list(df2['score'].astype(float))
	with open('int_{}.tsv'.format(cyt), 'w') as f:
		f.write('interactor\tscore\n')
		for i in range(len(list1)):
			f.write('{}\t{}\n'.format(list1[i], score1[i]))
		for i in range(len(list2)):
			f.write('{}\t{}\n'.format(list2[i], score2[i]))

with open('str_list.txt', 'r') as f:
	cyt_list = f.readline().split(',')

for cyt in cyt_list:
	extract_spec_interac(cyt)