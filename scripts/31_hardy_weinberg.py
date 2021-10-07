import pandas as pd

#reading snps data from tsv
df = pd.read_csv('./data/snp/snp_info.tsv', sep = '\t')
print(df)
#extracting frequencies of alleles as lists
p_list = list(df['freq_p'].astype(float))
q_list = list(df['freq_q'].astype(float))

#creating lists of frequencies for pp, qq, pq
qq = []
pp = []
pq = []

for p, q in zip(p_list, q_list):
	pp.append(p*p)
	qq.append(q*q)
	pq.append(2*p*q)

df['pp'] = pp
df['pq'] = pq
df['qq'] = qq

print(df)
df.to_csv('./data/snp/snp_hw.tsv', sep = '\t')