import requests
import pandas as pd
import urllib.error
import os

''' This module provides a class to work with the downloaded db
	of Gwas Catalog. I had to write this because I couldn't access
	(for whatever reason) the Gwas Catalog Rest API documentation.
	The class is initialitiated by downloading the GWAS Catalog
	and then extracting the various fields to create a Gwas Catalog
	object. Getter methods will be defined to obtain required information '''

url = 'https://www.ebi.ac.uk/gwas/api/search/downloads/alternative'
#catalog class
class GwasCatalog:
	def __init__(self, catalog):
		self.catalog = catalog

	#prints the catalog
	def showCatalog(self):
		print(self.catalog)

	#shows catalog columns indexes
	def showAttributes(self):
		for column in list(self.catalog.columns):
			print(column)

	#extracts single column from catalog as pandas Series object
	def getColumn(self, index):
		series = self.catalog[index]
		return series

	#search a column and returns all matching values rows
	def batchSearch(self, column_name, ids):
		results = self.catalog[self.catalog[column_name].isin(ids)]
		return results

	#same as batchSearch() but returns only selected columns (features)
	def batchRetrieve(self, column_name, ids, features):
		df = self.catalog[self.catalog[column_name].isin(ids)]
		series = [df[column_name]]
		for feature in features:
			series.append(df[feature])
		print(len(series))
		dataf = pd.DataFrame(series).T
		return dataf

#get catalog from local file
def getCatalog(path):
	df = pd.read_csv(path, sep = '\t', header = 0, low_memory = False)
	catalog = GwasCatalog(df)
	return catalog

#download latest catalog from EBI GwasCatalog
def updatedCatalog(filename, url = url, remove = True):
	try:
		print('Downloading GWAS catalog')
		print('Depending on your connection speed this may take up to some minutes...')
		f = requests.get(url).text
		print('Download Completed\n')
		with open(filename, 'w') as file:
			file.write(f)
		df = pd.read_csv(filename, sep = '\t', low_memory = False)
		if remove:
			os.remove(filename)
	except urllib.error.HTTPError:
		print('The URL is no longer valid or the server is unreachable.\n')
		url = input('Please insert a new URL: ')
		try:
			print('Downloading GWAS catalog')
			print('Depending on your connection speed this may take up to some minutes...')
			f = requests.get(url).text
			print('Download Completed\n')
			with open(filename, 'w') as file:
				file.write(f)
			df = pd.read_csv(filename, sep = '\t', low_memory = False)
			os.remove('catalog.tsv')
		except urllib.error.HTTPError:
			print('Operation couldn\'t be performed. Exiting...')
			exit()
	catalog = GwasCatalog(df)
	return catalog