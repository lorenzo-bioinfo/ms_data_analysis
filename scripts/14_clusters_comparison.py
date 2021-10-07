import pandas as pd

def getGroups(file, stops):
	ordered_ids = []
	with open(file) as f:
		for line in f:
			ordered_ids.append(int(line.strip()))
	mask = []
	for i in range(len(ordered_ids)):
		if ordered_ids[i] in stops:
			mask.append(-1)
		else:
			mask.append(0)
	groups = []
	gp = []
	for i in range(len(mask)):
		if mask[i] == 0:
			gp.append(ordered_ids[i])
		else:
			groups.append(gp)
			gp = [ordered_ids[i]]
	groups.append(gp)
	return groups

int_stops = [1672, 1474, 484, 1577, 1088, 698]
int_groups = getGroups('./data/cluster_groups/ordered_rr_ward.txt', int_stops)
print('int groups: ', len(int_groups))
norm_stops = [548, 1394, 1121, 1693, 1107]
norm_groups = getGroups('./data/cluster_groups/ordered_rr_norm.txt', norm_stops)
print('norm groups: ', len(norm_groups))
matrix = []
for a in int_groups:
	comp = []
	for b in norm_groups:
		comp.append(len(set(a) & set(b)) / len(a))
	matrix.append(comp)

rows = 'G1I,G2I,G3I,G4I,G5I,G6I,G7I'.split(',')
columns = 'G1N,G2N,G3N,G4N,G5N,G6N'.split(',')
mat = pd.DataFrame(matrix)
mat.index = rows
mat.columns = columns
print(mat)