import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from scipy.cluster import hierarchy

#getting cytokynes group for each cluster

colmap = ['darkgrey', 'darkgreen', 'navy']

clusters = []
for i in range(0, 7):
	with open('./data/cluster_groups/cyt_groups{}.txt'.format(i), 'r') as f:
		cluster = []
		group = []
		for line in f:
			clean_line = line.strip()
			if clean_line == '+':
				cluster.append(group)
				group = []
			else:
				group.append(clean_line)
		cluster.pop(0)
		clusters.append(cluster)

cyt_list = 'IL1B,IL2,IL4,IL5,IL6,IL7,CXCL8,IL10,IL12B,IL13,IL17A,CSF3,CSF2,IFNG,CCL2,CCL4,TNF,IL1RN,IL9,IL15,CCL11,FGF2,CXCL10,PDGFB,CCL5,VEGFA,CCL3'.split(',')
cyt_ord = ['IL1B', 'IL2', 'IL4', 'IL6', 'IL7', 'IL12B', 'IL17A', 'CSF2', 'IL15', 'FGF2', 'PDGFB', 'VEGFA', 'CCL3', 'CSF3', 'IL5', 'IL1RN', 'CXCL8', 'IL10', 'IL13', 'CCL11', 'CXCL10', 'CCL5', 'IL9', 'TNF', 'CCL4', 'CCL2', 'IFNG']
#cyt_ord = ['IL1B', 'IL2', 'IL4', 'IL6', 'IL7', 'IL5', 'IL1RN', 'CXCL8', 'IL10', 'IL13', 'IL12B', 'IL17A', 'CSF2', 'IL15', 'FGF2', 'PDGFB', 'VEGFA', 'CCL3', 'CSF3', 'CCL11', 'CXCL10', 'CCL5', 'IL9', 'TNF', 'CCL4', 'CCL2', 'IFNG']

#creating a 27x27 matrix which counts how many times each cytokine has clustered
#with all the other ones

zero_matrix = []
line = [0] * 27
for i in range(27):
	zero_matrix.append(line)

matrix = pd.DataFrame(zero_matrix, index = cyt_list, columns = cyt_list)

for cluster in clusters:
	for group in cluster:
		for el in group:
			for cyt in cyt_list:
				if cyt in group:
					matrix[el][cyt] += 1

#exporting matrix to file

matrix.to_csv('./data/cluster_groups/occ_matrix.tsv', sep = '\t')

#representing matrix in heatmap

sns.heatmap(matrix, cmap = 'mako')
#plt.savefig('./data/cluster_groups/heatmaps/occ_hm.png', dpi = 300)
plt.clf()

#generating a new matrix with cytokines in different order
#and using it to create a new heatmap

matrix_ord = matrix.reindex(index = cyt_ord, columns = cyt_ord)
sns.heatmap(matrix_ord, cmap = 'mako')
#plt.savefig('./data/cluster_groups/heatmaps/occ_hm_ord.png', dpi = 300)
plt.clf()

#Using clustering to take a look a the groups it forms

cluster_col = hierarchy.linkage(matrix.T, method="ward", metric="euclidean")
cluster_row = hierarchy.linkage(matrix, method="ward", metric="euclidean")
clusterfig = sns.clustermap(matrix, row_linkage = cluster_row, col_linkage = cluster_col)
index_col = clusterfig.dendrogram_col.reordered_ind
index_row = clusterfig.dendrogram_row.reordered_ind
plt.title('Cyt Clustering')
plt.savefig('./data/cluster_groups/heatmaps/occ_cluster.png', dpi = 300)
plt.clf()

#creating a 27x27 matrix which contains the mean distance of each
#cytokine from the other ones, in the ordered list generated by
#each cluster

zero_matrix = []
for i in range(27):
	zero_matrix.append(line)

matrix = pd.DataFrame(zero_matrix, index = cyt_list, columns = cyt_list)

#generating ordered lists from groups

ordered_lists = []
for cluster in clusters:
	ordered_list = []
	for group in cluster:
		ordered_list.extend(group)
	ordered_lists.append(ordered_list)

for ol in ordered_lists:
	for cyt_a in cyt_list:
		for cyt_b in cyt_list:
			dist = abs(ol.index(cyt_b) - ol.index(cyt_a))
			matrix[cyt_a][cyt_b] += dist

#dividing all cumulative distances by number of clusters
#to get mean distance

def dividi(x):
	return(x/7)

matrix_ok = matrix.applymap(dividi)
matrix_ok.to_csv('./data/cluster_groups/dist_matrix.tsv', sep = '\t')
sns.heatmap(matrix_ok, cmap = 'mako')
#plt.savefig('./data/cluster_groups/heatmaps/dist_hm.png', dpi = 300)
plt.clf()

matrix_ord = matrix.reindex(index = cyt_ord, columns = cyt_ord)
sns.heatmap(matrix_ord, cmap = 'mako')
#plt.savefig('./data/cluster_groups/heatmaps/dist_hm_ord.png', dpi = 300)
plt.clf()

cluster_col = hierarchy.linkage(matrix_ok.T, method="ward", metric="euclidean")
cluster_row = hierarchy.linkage(matrix_ok, method="ward", metric="euclidean")
clusterfig = sns.clustermap(matrix_ok, row_linkage = cluster_row, col_linkage = cluster_col)
index_col = clusterfig.dendrogram_col.reordered_ind
index_row = clusterfig.dendrogram_row.reordered_ind
plt.title('Cyt Clustering')
plt.savefig('./data/cluster_groups/heatmaps/dist_cluster.png', dpi = 300)
plt.clf()