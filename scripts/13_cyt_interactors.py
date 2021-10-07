#In this script I will use the data extracted from STRING
#(check ./data/string_networks) to check how many interactors
#have been reported for each cytokine, and how many interactors
#are common to ALL of the cytokines

#importing libraries
import pandas as pd

#Getting list of cytokines

cyt_list = 'IL1B,IL2,IL4,IL5,IL6,IL7,CXCL8,IL10,IL12B,IL13,IL17A,CSF3,CSF2,IFNG,CCL2,CCL4,TNF,IL1RN,IL9,IL15,CCL11,FGF2,CXCL10,PDGFB,CCL5,VEGFA,CCL3'.split(',')

#Reading files generated with strings, in order to obtain a dictionary
#in the form '{cyt name = [list of interactors]}'
#In the mean time, I will create a list containing ALL of
#the interactors, repeated only once

#initializing dictionary
cyt_interactors = {}
#initializing interactors list
interactors = []
#extracting interactors from files
for cyt in cyt_list:
	df = pd.read_csv('./data/string_networks/int_{}.tsv'.format(cyt), sep = '\t')
	int_list = list(df['interactor'])
	#adding entry to dictionary
	#every interactor is repeated twice for each
	#cytokine, so i'll cast the list to set and
	#then to list again, in order to remove duplicates
	cyt_interactors[cyt] = list(set(int_list))
	#adding interactors to list
	interactors.extend(int_list)


#getting number of interactors for each cytokine
for cyt in cyt_list:
	print(cyt, len(cyt_interactors[cyt]))
#about 490 for each cytokine, with a drop at 213 for IL12B

#cleaning interactors list from duplicates
interactors_clean = list(set(interactors))
print('Total number of interactors: {}'.format(len(interactors_clean))) #2726
#exporting interactors into a file
with open('./data/string_networks/interactors_list.txt', 'w') as f:
	for interactor in interactors_clean:
		f.write(interactor + '\n')

#Looking for an efficient way to get only the common interactors

common_cyt = set(cyt_interactors['IL1B']) #I have to start somewhere!
for cyt in cyt_list:
	cyt_int = set(cyt_interactors[cyt])
	common_cyt = common_cyt & cyt_int
print('Common interactors: {}'.format(len(common_cyt)))

#And now getting a matrix which shows how many interactors
#each cytokine has in common with all the other ones
com_matrix = []
for cyt_a in cyt_list:
	com = []
	int_lst = set(cyt_interactors[cyt_a])
	for cyt_b in cyt_list:
		com.append(len(int_lst & set(cyt_interactors[cyt_b])))
	com_matrix.append(com)

#turning the matrix into a dataframe and exporting it to a tsv file 

com_df = pd.DataFrame(com_matrix)
com_df.index = cyt_list
com_df.columns = cyt_list
print(com_df)
com_df.to_csv('./data/string_networks/common_int_matrix.tsv', sep = '\t')

#with the method above I'm getting only 6 (SIX) common interactors.
#Trying to get rid of IL12B to see if it gets better

cyt_list.pop(cyt_list.index('IL12B'))
common_cyt = set(cyt_interactors['IL1B']) #I have to start somewhere!
for cyt in cyt_list:
	cyt_int = set(cyt_interactors[cyt])
	common_cyt = common_cyt & cyt_int
print('Common interactors without IL12B: {}'.format(len(common_cyt)))

#22 now, better but probably not enough :(

#All the interactors have been fed to David for functional annotation.
#A table has been downloaded containing the protein name and all
#the GO:BP terms associated with it.
#Trying to read that file (I can see problems comin :'( )

david_df = pd.read_csv('./data/string_networks/david_GO_BP.txt', sep = '\t')
print(david_df.describe())
print(david_df)

#looks like this worked. Now I'll try to get a dictionary in the form:
# {interactor_name: [list of associated GO:BP terms]}
#At the same time I'll get a list of all the terms

#this is used to iterate the david_df rows
index = list(david_df['ID'])

go_bp_dict = {}
go_terms = []

for identifier in index:
	go = list(david_df[david_df['ID'] == identifier]['GOTERM_BP_DIRECT'].str.split(','))[0]
	go_bp_dict[identifier] = go
	for string in go:
		try:
			go_terms.append(string.split('~')[1].strip())
		except IndexError:
			pass
#cleaning the terms list
go_terms_clean = list(set(go_terms))
#writing the list to a file
with open('./data/string_networks/go_terms.txt', 'w') as f:
	for entry in go_terms_clean:
		f.write(entry + '\n')

#checking results:
print('Number of GO terms: {}'.format(len(go_terms_clean)))

#exporting interactors AND their annotation in a file

with open('./data/string_networks/annot_interactors.txt', 'w') as f:
	for interactor in go_bp_dict:
		f.write('{}^'.format(interactor))
		f.write('$'.join(go_bp_dict[interactor]) + '\n')