import GwasCatalogDigger as gcd
import pandas as pd

catalog = gcd.getCatalog('./data/snp/gwas_catalog.tsv')
print(catalog)
print(type(catalog))
catalog.showCatalog()
catalog.showAttributes()
snp_col = catalog.getColumn('SNPS')
print(snp_col)

clinvar = pd.read_csv('data/snp/snp_clinvar.tsv', sep = '\t')
snp_list = list(clinvar['snp_id'])

features = ['DISEASE/TRAIT', 'CHR_ID', 'STRONGEST SNP-RISK ALLELE', 'P-VALUE', 'PUBMEDID']

dataf = catalog.batchRetrieve('SNPS', snp_list, features)
new_col = []
new_col2 = []
new_col3 = []
new_col4 = []
snp_list = list(dataf['SNPS'])
for snp in snp_list:
    new_col.append(clinvar[clinvar['snp_id'] == snp]['gene_consequence1'].values[0])
    new_col2.append(clinvar[clinvar['snp_id'] == snp]['gene_consequence2'].values[0])
    new_col3.append(clinvar[clinvar['snp_id'] == snp]['gene_1'].values[0])
    new_col4.append(clinvar[clinvar['snp_id'] == snp]['gene_2'].values[0])

dataf['NCBI_gene1'] = new_col3
dataf['NCBI_info1'] = new_col
dataf['NCBI_gene2'] = new_col4
dataf['NCBI_info2'] = new_col2
print(dataf)
dataf.to_csv('./data/snp/gwas_info.tsv', sep = '\t')