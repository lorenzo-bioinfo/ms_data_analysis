import pandas as pd

#using p-values generated in ks_testing to
#apply a benjamini-hochberg correction

def benjamini_hochberg(dataframe, dataframe_column):
	#dataframe_column is the name or index of the col containing p-values
	dataframe.sort_values(by = dataframe_column, inplace = True)
	p_values = list(dataframe[dataframe_column].astype(float))
	q_values = []
	n = len(p_values)
	for i, p_value in enumerate(p_values):
		q_values.append((p_value * n)/(i+1))
	#adding q_values column to dataframe
	dataframe['q_values'] = q_values
	return dataframe

df_PP = pd.read_csv('./consistency_ks_cyt/ks_PP.tsv', sep = '\t', header = 0, index_col = 0)
df_SP = pd.read_csv('./consistency_ks_cyt/ks_SP.tsv', sep = '\t', header = 0, index_col = 0)
df_RR = pd.read_csv('./consistency_ks_cyt/ks_RR.tsv', sep = '\t', header = 0, index_col = 0)

df_PP_q = benjamini_hochberg(df_PP, 'p_value')
df_PP_q.to_csv('./consistency_ks_cyt/PP_bh_corr.tsv', sep = '\t')
df_SP_q = benjamini_hochberg(df_SP, 'p_value')
df_SP_q.to_csv('./consistency_ks_cyt/SP_bh_corr.tsv', sep = '\t')
df_RR_q = benjamini_hochberg(df_RR, 'p_value')
df_RR_q.to_csv('./consistency_ks_cyt/RR_bh_corr.tsv', sep = '\t')