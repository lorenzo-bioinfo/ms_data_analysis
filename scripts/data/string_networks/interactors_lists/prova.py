with open('full_interactors.txt', 'r') as f:
	lista = []
	for line in f:
		lista.append(line.strip())
print(len(lista), len(set(lista)))