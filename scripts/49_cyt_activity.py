import pandas as pd
from scipy import stats
from numpy import median

cyt_list = 'IL1B,IL2,IL4,IL5,IL6,IL7,CXCL8,IL10,IL12B,IL13,IL17A,CSF3,CSF2,IFNG,CCL2,CCL4,TNF,IL1RN,IL9,IL15,CCL11,FGF2,CXCL10,PDGFB,CCL5,VEGFA,CCL3'.split(',')
df = pd.read_excel('../database/db.xlsx', sheet_name = 'SM NM', usecols = 'A,F:AF,CB,CN,DD')
df.columns = 'patient_id,IL1B,IL2,IL4,IL5,IL6,IL7,CXCL8,IL10,IL12B,IL13,IL17A,CSF3,CSF2,IFNG,CCL2,CCL4,TNF,IL1RN,IL9,IL15,CCL11,FGF2,CXCL10,PDGFB,CCL5,VEGFA,CCL3,class_int,act,cort_therapy'.split(',')
df_rr = df[df['class_int'] == 3]
df_rr_filled = df_rr.fillna({'cort_therapy' : 1, 'act' : 0})
df_ok = df_rr_filled[df_rr_filled['cort_therapy'] == 1].dropna()

df_noact = df_ok[(df_ok['act'] == 0) | (df_ok['act'] == 'NO') | (df_ok['act'] == 'NA')]
df_act = df_ok[(df_ok['act'] == 1) | (df_ok['act'] == 2) | (df_ok['act'] == 3)]

with open('./data/activity/cyt_act.tsv', 'w') as f:
	f.write('cyt\tact_med\tnoact_med\tks_stat\tp_val\n')
	for cyt in cyt_list:
		med_noact = df_noact[cyt].astype(float).median()
		med_act = df_act[cyt].astype(float).median()
		s1 = list(df_noact[cyt].astype(float))
		s2 = list(df_act[cyt].astype(float))
		stat = stats.ks_2samp(s1, s2)
		f.write(f'{cyt}\t{med_act}\t{med_noact}\t{stat[0]}\t{stat[1]}\n')

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
	dataframe.sort_values('q_values', inplace = True)
	return dataframe

df_stat = pd.read_csv('./data/activity/cyt_act.tsv', sep = '\t')
print(df_stat)
df_bh = benjamini_hochberg(df_stat, 'p_val')
df_bh.to_csv('./data/activity/cyt_act_bh.tsv', sep = '\t')