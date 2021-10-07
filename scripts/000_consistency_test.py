import pandas as pd
from scipy import stats

cyt_list = 'IL-1b,IL-2,IL-4,IL-5,IL-6,IL-7,IL-8,IL-10,IL-12,IL-13,IL-17,G-CSF,GM-CSF,IFN-g,MCP-1,MIP-1b,TNF-a,IL-1ra,IL-9,IL-15,Eoxatin,FGC_basic,IP-10,PDGF-bb,RANTES,VEGF,MIP-1a'.split(',')
old_ct = []
with open('./data/old_ctrl.txt', 'r') as f:
	for line in f:
		old_ct.append(int(line.strip()))

old_rr = []
with open('./data/old_rr.txt') as f:
	for line in f:
		old_rr.append(int(line.strip()))

#getting df from csv with pre-PL therapy infos

df = pd.read_excel('../database/db.xlsx', sheet_name = 'SM NM', usecols = 'A,F:AF,CB,DD')
df.columns = ('patient_id,IL-1b,IL-2,IL-4,IL-5,IL-6,IL-7,IL-8,IL-10,IL-12,IL-13,IL-17,G-CSF,GM-CSF,IFN-g,MCP-1,MIP-1b,TNF-a,IL-1ra,IL-9,IL-15,Eoxatin,FGC_basic,IP-10,PDGF-bb,RANTES,VEGF,MIP-1a,class_int,cort_therapy').split(',')

#This lines were useless, keeping them here fo precaution's sake
#replacing NaN values with 1 (no therapy)
#df['cort_therapy'].fillna(1, inplace = True) #this doesn't work!!!

#creating sub-dfs for each class (no therapy before pl)
df_ctrl = df[df['class_int'] == 6]
df_ctrl_filled = df_ctrl.fillna({'cort_therapy' : 1})
df_ctrl_full = df_ctrl_filled[df_ctrl_filled['cort_therapy'] == 1].dropna()
df_ctrl_noc = df_ctrl_full[df_ctrl_full['patient_id'].isin(old_ct)]
df_pp = df[df['class_int'] == 5]
df_pp_filled = df_pp.fillna({'cort_therapy' : 1})
df_pp_noc = df_pp_filled[df_pp_filled['cort_therapy'] == 1].dropna()
df_sp = df[df['class_int'] == 4]
df_sp_filled = df_sp.fillna({'cort_therapy' : 1})
df_sp_noc = df_sp_filled[df_sp_filled['cort_therapy'] == 1].dropna()
df_rr = df[df['class_int'] == 3]
df_rr_filled = df_rr.fillna({'cort_therapy' : 1})
df_rr_full = df_rr_filled[df_rr_filled['cort_therapy'] == 1].dropna()
df_rr_noc = df_rr_full[df_rr_full['patient_id'].isin(old_rr)]
print(len(df_ctrl_noc), len(df_pp_noc), len(df_sp_noc), len(df_rr_noc))

missing_ct = [x for x in old_ct if x not in list(df_ctrl_noc['patient_id'])]
missing_rr = [x for x in old_rr if x not in list(df_rr_noc['patient_id'])]
print('Missing ctrl:')
print(missing_ct, len(missing_ct))
print('Missing RR:')
print(missing_rr, len(missing_rr))

#using Kolmogorov-Smirnov 2 sample test to see if there are significative
#differences between cytokines levels in control and patients

datafs = [('PP', df_pp_noc), ('SP', df_sp_noc), ('RR', df_rr_noc)]

for dataf in datafs:
	with open('./consistency_ks_cyt/ks_{}.tsv'.format(dataf[0]), 'w') as f:
		f.write('cyt\tks_stat\tp_value\tctrl_median\tpatient_median\n')
		for cyt in cyt_list:
			f.write(cyt + '\t')
			ctrl_median = df_ctrl_noc[cyt].astype(float).median()
			class_median = dataf[1][cyt].astype(float).median()
			stat = stats.ks_2samp(df_ctrl_noc[cyt].astype(float), dataf[1][cyt].astype(float))
			f.write('{:0.6f}\t'.format(stat[0]))
			f.write('{:0.6f}\t'.format(stat[1]))
			f.write('{:0.6f}\t'.format(ctrl_median))
			f.write('{:0.6f}\n'.format(class_median))

#results were stored in data/ks_test
#Interpretation:
	#ks_stat: maximum absolute distance found between samples
	#p_value: just a regular p_value
	#the closest ks_stat is to 0, the more the two samples are likely to belong
	#to the same distribution

#cleaning data from outliers before repeating the test:
	#criteria: values that differ from mean more than 2 std are considered outliers

for dataf in datafs:
	with open('./consistency_ks_cyt/ks_clean_{}.tsv'.format(dataf[0]), 'w') as f:
		f.write('cyt\tks_stat\tp_value\tctrl_median\tpatient_median\n')
		for cyt in cyt_list:
			f.write(cyt + '\t')
			series_ctrl = df_ctrl_noc[cyt].astype(float)
			ctrl_mean = series_ctrl.mean() 
			ctrl_std = series_ctrl.std()
			series_ctrl.drop(series_ctrl.index[(series_ctrl <= ctrl_mean - ctrl_std * 2) | (series_ctrl >= ctrl_mean + ctrl_std * 2)], inplace = True)
			ctrl_median = series_ctrl.median()
			series_class = dataf[1][cyt].astype(float)
			class_mean = series_class.mean()
			class_std = series_class.std()
			series_class.drop(series_class.index[(series_class <= class_mean - class_std * 2) | (series_class >= class_mean + class_std * 2)], inplace = True)
			class_median = series_class.median()
			stat = stats.ks_2samp(series_ctrl, series_class)
			f.write('{:0.6f}\t'.format(stat[0]))
			f.write('{:0.6f}\t'.format(stat[1]))
			f.write('{:0.6f}\t'.format(ctrl_median))
			f.write('{:0.6f}\n'.format(class_median))

#obtain snips from db (url, biopython, chose)
#use scipy.cluster.hierarchy to divide in groups and see stuff
print('################### CITOCHINE ERRERRE ###################')
for cyt in cyt_list:
	print(cyt, '---', len(list(df_rr_noc[cyt])))

print('################### CITOCHINE CONTROLLI ###################')
for cyt in cyt_list:
	print(cyt, '---', len(list(df_ctrl_noc[cyt])))

missing_df = df[df['patient_id'].isin(missing_rr)]
print(missing_df)
for el in missing_rr:
	print(missing_df[missing_df['patient_id'] == el])