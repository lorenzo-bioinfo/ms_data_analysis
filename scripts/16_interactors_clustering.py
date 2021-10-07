import pandas as pd
from scipy.cluster import hierarchy
from matplotlib import pyplot as plt
import seaborn as sns

#reading list of interactors and their GO:BP annotations

interactors_dict = {}

with open('./data/string_networks/annot_interactors.txt', 'r') as f:
	for line in f:
		interactor, terms_str = line.strip().split('^')
		terms_list = terms_str.split('$')
		terms_clean = []
		for term in terms_list:
			try:
				terms_clean.append(term.split('~')[1].strip())
			except IndexError:
				pass
		interactors_dict[interactor] = terms_clean

#reading the list of go terms "associated" with multiple sclerosis

with open('./data/string_networks/go_terms_clean.txt', 'r') as f:
	ms_terms = []
	for line in f:
		ms_terms.append(line.strip())

#reading the list of go terms "associated" with myelination

with open('./data/string_networks/go_terms_myelin.txt', 'r') as f:
	myelin_terms = []
	for line in f:
		myelin_terms.append(line.strip())

#selecting only interactors which have at least one term
#also present in ms_terms

ms_interactors = [] #1658 interactors

f = open('./data/string_networks/interactors_lists/ms_interactors.txt', 'w')
for interac in interactors_dict:
	for term in interactors_dict[interac]:
		if term in ms_terms:
			ms_interactors.append(interac)
			f.write('{}\n'.format(interac))
			break
f.close()

#selecting only interactors which have at least one term
#also present in myelin_terms

myelin_interactors = []

f = open('./data/string_networks/interactors_lists/myelin_interactors.txt', 'w')
for interac in interactors_dict:
	for term in interactors_dict[interac]:
		if term in myelin_terms:
			myelin_interactors.append(interac)
			f.write('{}\n'.format(interac))
			break
f.close()

#importing table of cytokines and their interactors

cyt_list = 'IL1B,IL2,IL4,IL5,IL6,IL7,CXCL8,IL10,IL12B,IL13,IL17A,CSF3,CSF2,IFNG,CCL2,CCL4,TNF,IL1RN,IL9,IL15,CCL11,FGF2,CXCL10,PDGFB,CCL5,VEGFA,CCL3'.split(',')
cyt_interactors = {}
for cyt in cyt_list:
	df = pd.read_csv('./data/string_networks/int_{}.tsv'.format(cyt), sep = '\t')
	int_list = list(df['interactor'])
	scores_list = list(df['score'].astype(float))
	length = range(len(int_list))
	int_scores = []
	for i in length:
		int_scores.append((int_list[i], scores_list[i]))
	cyt_interactors[cyt] = list(set(int_scores))

#creating the matrix for ms_interactors that will be used to generate clusters
#for each cytokine, the score of the interaction with every interactor
#in ms_interactors will be reported. In case of no interaction, -1 will
#be used

int_matrix = []

for cyt in cyt_list:
	matrix_row = []
	cyt_intscores = cyt_interactors[cyt]
	cyt_int = []
	cyt_score = []
	for interac, score in cyt_intscores:
		cyt_int.append(interac)
		cyt_score.append(score)
	for ms_int in ms_interactors:
		if ms_int in cyt_int:
			matrix_row.append(cyt_score[cyt_int.index(ms_int)])
		else:
			matrix_row.append(-1)
	int_matrix.append(matrix_row)

int_df = pd.DataFrame(int_matrix)
int_df.index = cyt_list
int_df.columns = ms_interactors

#creating the matrix for myelin_interactors that will be used to generate clusters
#for each cytokine, the score of the interaction with every interactor
#in myelin_interactors will be reported. In case of no interaction, -1 will
#be used

myelin_matrix = []

for cyt in cyt_list:
	matrix_row = []
	cyt_intscores = cyt_interactors[cyt]
	cyt_int = []
	cyt_score = []
	for interac, score in cyt_intscores:
		cyt_int.append(interac)
		cyt_score.append(score)
	for ms_int in myelin_interactors:
		if ms_int in cyt_int:
			matrix_row.append(cyt_score[cyt_int.index(ms_int)])
		else:
			matrix_row.append(-1)
	myelin_matrix.append(matrix_row)

myelin_df = pd.DataFrame(myelin_matrix)
myelin_df.index = cyt_list
myelin_df.columns = myelin_interactors

#generating cytokines clusters using int_df

cluster_col = hierarchy.linkage(int_df.T, method="ward", metric="euclidean")
cluster_row = hierarchy.linkage(int_df, method="ward", metric="euclidean")
clusterfig = sns.clustermap(int_df, row_linkage = cluster_row, col_linkage = cluster_col, yticklabels = True, figsize = (10, len(int_df)/4), cmap = 'mako')
index_col = clusterfig.dendrogram_col.reordered_ind #cytokines
index_row = clusterfig.dendrogram_row.reordered_ind #patients
plt.savefig('../plots/cyt_interactors_clusters.png', dpi = 300)
plt.clf()

#generating cytokines clusters using myelin_df

cluster_col = hierarchy.linkage(myelin_df.T, method="ward", metric="euclidean")
cluster_row = hierarchy.linkage(myelin_df, method="ward", metric="euclidean")
clusterfig = sns.clustermap(myelin_df, row_linkage = cluster_row, col_linkage = cluster_col, yticklabels = True, figsize = (10, len(myelin_df)/4), cmap = 'mako')
index_col = clusterfig.dendrogram_col.reordered_ind #cytokines
index_row = clusterfig.dendrogram_row.reordered_ind #patients
plt.savefig('../plots/cyt_interactors_myelination.png', dpi = 300)
plt.clf()