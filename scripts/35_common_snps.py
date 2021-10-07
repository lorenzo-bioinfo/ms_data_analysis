#Isolating snps which genotypes differ both in RR vs Euro and
#RR vs controls. This SHOULD give me RR 'specific' snps

import pandas as pd

df_euro = pd.read_csv('./data/snp/tests/rr_vs_euro_bh.tsv', sep = '\t', header = 0, index_col = 0)
df_ctrl = pd.read_csv('./data/snp/tests/ctrl_vs_rr_bh.tsv', sep = '\t', header = 0, index_col = 0)
print(df_euro, df_ctrl)

df_euro_sign = df_euro[df_euro['q_values'] < 0.05]
df_ctrl_sign = df_ctrl[df_ctrl['q_values'] < 0.05]

euro_set = set(list(df_euro_sign['snp_id']))
ctrl_set = set(list(df_ctrl_sign['snp_id']))

common_snps = euro_set & ctrl_set
print(f'Euro --> {euro_set} | {len(euro_set)}')
print(f'Ctrl --> {ctrl_set} | {len(ctrl_set)}')
print(f'Common --> {common_snps} |{len(common_snps)}')