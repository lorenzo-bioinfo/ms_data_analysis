import pandas as pd
from scipy import stats

cyt_list = 'IL1B,IL2,IL4,IL5,IL6,IL7,CXCL8,IL10,IL12B,IL13,IL17A,CSF3,CSF2,IFNG,CCL2,CCL4,TNF,IL1RN,IL9,IL15,CCL11,FGF2,CXCL10,PDGFB,CCL5,VEGFA,CCL3'.split(',')

#getting df from csv with pre-PL therapy infos

df = pd.read_excel('../database/db.xlsx', sheet_name = 'SM NM', usecols = 'A,F:AF,CB,DD')
df.columns = ('patient_id,IL1B,IL2,IL4,IL5,IL6,IL7,CXCL8,IL10,IL12B,IL13,IL17A,CSF3,CSF2,IFNG,CCL2,CCL4,TNF,IL1RN,IL9,IL15,CCL11,FGF2,CXCL10,PDGFB,CCL5,VEGFA,CCL3,class_int,cort_therapy').split(',')

#This lines were useless, keeping them here fo precaution's sake
#replacing NaN values with 1 (no therapy)
#df['cort_therapy'].fillna(1, inplace = True) #this doesn't work!!!

#creating sub-dfs for each class (no therapy before pl)
df_ctrl = df[df['class_int'] == 6]
df_ctrl_filled = df_ctrl.fillna({'cort_therapy' : 1})
df_ctrl_noc = df_ctrl_filled[df_ctrl_filled['cort_therapy'] == 1]
df_pp = df[df['class_int'] == 5]
df_pp_filled = df_pp.fillna({'cort_therapy' : 1})
df_pp_noc = df_pp_filled[df_pp_filled['cort_therapy'] == 1]
df_sp = df[df['class_int'] == 4]
df_sp_filled = df_sp.fillna({'cort_therapy' : 1})
df_sp_noc = df_sp_filled[df_sp_filled['cort_therapy'] == 1]
df_rr = df[df['class_int'] == 3]
df_rr_filled = df_rr.fillna({'cort_therapy' : 1})
df_rr_noc = df_rr_filled[df_rr_filled['cort_therapy'] == 1]

print(len(df_ctrl_noc), len(df_pp_noc), len(df_sp_noc), len(df_rr_noc))

#exporting correct df's to avoid further problems

df_ctrl_noc.to_csv('./data/dataframes/ctrl_noc.csv')
df_pp_noc.to_csv('./data/dataframes/pp_noc.csv')
df_sp_noc.to_csv('./data/dataframes/sp_noc.csv')
df_rr_noc.to_csv('./data/dataframes/rr_noc.csv')

#using Kolmogorov-Smirnov 2 sample test to see if there are significative
#differences between cytokines levels in control and patients

datafs = [('PP', df_pp_noc), ('SP', df_sp_noc), ('RR', df_rr_noc)]

for dataf in datafs:
	with open('./data/ks_test/ks_{}.tsv'.format(dataf[0]), 'w') as f:
		f.write('cyt\tks_stat\tp_value\tctrl_median\tpatient_median\tN_ct\tN_pat\n')
		for cyt in cyt_list:
			f.write(cyt + '\t')
			ctrl_median = df_ctrl_noc[cyt].dropna().astype(float).median()
			class_median = dataf[1][cyt].dropna().astype(float).median()
			stat = stats.ks_2samp(df_ctrl_noc[cyt].dropna().astype(float), dataf[1][cyt].dropna().astype(float))
			n_ct = len(df_ctrl_noc[cyt].dropna())
			n_pat = len(dataf[1][cyt].dropna())
			f.write('{:0.6f}\t'.format(stat[0]))
			f.write('{:0.6f}\t'.format(stat[1]))
			f.write('{:0.6f}\t'.format(ctrl_median))
			f.write('{:0.6f}\t'.format(class_median))
			f.write('{}\t'.format(n_ct))
			f.write('{}\n'.format(n_pat))

#results were stored in data/ks_test
#Interpretation:
	#ks_stat: maximum absolute distance found between samples
	#p_value: just a regular p_value
	#the closest ks_stat is to 0, the more the two samples are likely to belong
	#to the same distribution

for dataf in datafs:
	with open('./data/ks_test/ks_clean_{}.tsv'.format(dataf[0]), 'w') as f:
		f.write('cyt\tks_stat\tp_value\tctrl_median\tpatient_median\tN_ct\tN_pat\n')
		for cyt in cyt_list:
			f.write(cyt + '\t')
			print('################# {} #################'.format(cyt))
			#dropping outliers from controls
			series_ctrl = df_ctrl_noc[cyt].dropna().astype(float)
			ctrl_lenout = (len(series_ctrl))
			norm_pval = stats.shapiro(series_ctrl)[1]
			if norm_pval > 0.05:
				ctrl_mean = series_ctrl.mean() 
				ctrl_std = series_ctrl.std()
				print('NORMAL CONTROLS!!!')
				series_ctrl.drop(series_ctrl.index[(series_ctrl < ctrl_mean - ctrl_std * 2) | (series_ctrl > ctrl_mean + ctrl_std * 2)], inplace = True)
			else:
				thresh = series_ctrl.quantile(0.95)
				series_ctrl.drop(series_ctrl.index[(series_ctrl > thresh)], inplace = True)
			ctrl_median = series_ctrl.median()
			ctrl_lenin = (len(series_ctrl))
			print('###CTRL:')
			print('Before: {}\t\t After: {}\n'.format(ctrl_lenout, ctrl_lenin))

			#dropping outliers from patients
			series_class = dataf[1][cyt].dropna().astype(float)
			class_lenout = len(series_class)
			norm_pval = stats.shapiro(series_class)[1]
			if norm_pval > 0.05:
				class_mean = series_class.mean()
				class_std = series_class.std()
				print('NORMAL PATIENTS!!!')
				series_class.drop(series_class.index[(series_class < class_mean - class_std * 2) | (series_class > class_mean + class_std * 2)], inplace = True)
			else:
				thresh = series_class.quantile(0.95)
				series_class.drop(series_class.index[series_class > thresh], inplace = True)
			class_median = series_class.median()
			class_lenin = len(series_class)
			print('###{}:'.format(dataf[0]))
			print('Before: {}\t\t After: {}\n'.format(class_lenout, class_lenin))

			#testing and exporting results
			stat = stats.ks_2samp(series_ctrl, series_class)
			f.write('{:0.8f}\t'.format(stat[0]))
			f.write('{:0.8f}\t'.format(stat[1]))
			f.write('{:0.8f}\t'.format(ctrl_median))
			f.write('{:0.8f}\t'.format(class_median))
			f.write('{}\t'.format(ctrl_lenin))
			f.write('{}\n'.format(class_lenin))
