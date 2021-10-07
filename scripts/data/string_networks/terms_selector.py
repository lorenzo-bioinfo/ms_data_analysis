#with this script I'm going to select only the GO terms
#which show a match with a word in the word list.
#In this way I hope to get ALL the terms relevant for
#multiple sclerosis, plus some clutter that will be easier to clean
#by hand

#getting word list
with open('./word_list.txt', 'r') as f:
	word_list = f.readline().strip().split(',')

#getting myelination word list
with open('./word_list_myelin.txt', 'r') as f:
	word_list_myelin = f.readline().strip().split(',')

#getting dictionary of interactors plus theri GO:BP annotation
int_dict = {}
with open('annot_interactors.txt', 'r') as f:
	for line in f:
		interactor, annot_str = line.strip().split('^')
		annot = annot_str.split('$')
		int_dict[interactor] = annot

#getting only terms that match the word list:
terms = []
myel_terms = []
for key in int_dict:
	for entry in int_dict[key]:
		try:
			annot = entry.split('~')[1].strip().split(' ')
			for word in annot:
				if word in word_list:
					terms.append(' '.join(annot))
				if word in word_list_myelin:
					myel_terms.append(' '.join(annot))
		except IndexError:
			for word in annot:
				if word in word_list:
					terms.append(' '.join(annot))
				if word in word_list_myelin:
					myel_terms.append(' '.join(annot))

#cleaning from duplicates:
terms_clean = list(set(terms))
myel_terms_clean = list(set(myel_terms))
print(terms_clean, len(terms_clean))
print(myel_terms_clean, len(myel_terms_clean))

i = 0

#exporting selected terms in a file
with open('./go_terms_clean.txt', 'w') as f:
	for term in terms_clean:
		f.write('{}\n'.format(term))
with open('./go_terms_myelin.txt', 'w') as f:
	for term in myel_terms_clean:
		f.write('{}\n'.format(term))