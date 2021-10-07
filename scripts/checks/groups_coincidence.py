import pandas as pd

with open('./../data/cluster_groups/ward_groups.txt', 'r') as f:
	wg = []
	for line in f:
		wg.append(line.strip().split(','))
print(len(wg))
with open('./../data/cluster_groups/interactors_groups.txt', 'r') as f:
	ig = []
	for line in f:
		ig.append(line.strip().split(','))
print(len(ig))
matrix = []
for wgg in wg:
	row = []
	for igg in ig:
		count = 0
		for el in igg:
			if el in wgg:
				count += 1
		count_ok = count/len(igg)
		row.append(count_ok)
	matrix.append(row)

df = pd.DataFrame(matrix)
df.index = [1, 2, 3, 4, 5, 6]
df.columns = [1, 2, 3, 4, 5]
print(df)