import pandas as pd

df = pd.read_csv('../data/snp/gwas_info.tsv', sep = '\t')
print(df)
snps = 'rs11913721,rs1474347,rs11568817,rs1554973,rs1799880,rs2069763,rs3757351,rs2227306,rs2069772,rs7301328,rs7598440,rs743409,rs4680,rs12722489,rs2227284,rs1800797,rs2283792,rs4251961,rs10767664,rs2069812,rs1927914,rs7131056'.split(',')
print(len(snps))
df_ok = df[df['SNPS'].isin(snps)]
print(df_ok)
df_ok.drop('P-VALUE', axis = 1, inplace = True)
df_ok.drop('STRONGEST SNP-RISK ALLELE', axis = 1, inplace = True)
df_ok.drop('PUBMEDID', axis = 1, inplace = True)
df_ok.to_csv('snps_ok_info_full.tsv', sep = '\t')