import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

cyt_list = 'IL1B,IL2,IL4,IL5,IL6,IL7,CXCL8,IL10,IL12B,IL13,IL17A,CSF3,CSF2,IFNG,CCL2,CCL4,TNF,IL1RN,IL9,IL15,CCL11,FGF2,CXCL10,PDGFB,CCL5,VEGFA,CCL3'.split(',')

#getting df from csv with pre-PL therapy infos

df = pd.read_excel('../database/db.xlsx', sheet_name = 'SM NM', usecols = 'A,F:AF,CB,DD')
df.columns = ('patient_id,IL1B,IL2,IL4,IL5,IL6,IL7,CXCL8,IL10,IL12B,IL13,IL17A,CSF3,CSF2,IFNG,CCL2,CCL4,TNF,IL1RN,IL9,IL15,CCL11,FGF2,CXCL10,PDGFB,CCL5,VEGFA,CCL3,class_int,cort_therapy').split(',')

#This lines were useless, keeping them here fo precaution's sake
#replacing NaN values with 1 (no therapy)
#df['cort_therapy'].fillna(1, inplace = True)

#creating sub-dfs for each class (no therapy before pl)
df_ctrl_noc = df[(df['class_int'] == 6) & (df['cort_therapy'] == 1)].dropna()
df_pp_noc = df[(df['class_int'] == 5) & (df['cort_therapy'] == 1.0)].dropna()
df_sp_noc = df[(df['class_int'] == 4) & (df['cort_therapy'] == 1.0)].dropna()
df_rr_noc = df[(df['class_int'] == 3) & (df['cort_therapy'] == 1.0)].dropna()

print(len(df_ctrl_noc), len(df_pp_noc), len(df_sp_noc), len(df_rr_noc))

#plotting distributions of each cytokine for each class against control

sns.distributions._has_statsmodels = False #needed to avoid kde error coming from sns using statsmodel

#CTRL VS PP
for cyt in cyt_list:
	plt.title('{} - PP vs CTRL(No CS Therapy)\nN = {}'.format(cyt, len(df_pp_noc)))
	sns.displot(df_ctrl_noc[cyt], color = 'grey')
	sns.displot(df_pp_noc[cyt], color = 'darkgreen')
	plt.legend(['Control', 'PP'])
	plt.xlabel('{} levels'.format(cyt))
	plt.savefig('./../plots/ctrl_pp/{}_noc.png'.format(cyt), dpi = 300)
	print('Saved ctrl_pp/{}'.format(cyt))
	plt.clf()

#CTRL VS SP
for cyt in cyt_list:
	plt.title('{} - SP vs CTRL(No CS Therapy\nN = {})'.format(cyt, len(df_sp_noc)))
	sns.displot(df_ctrl_noc[cyt], color = 'grey')
	sns.displot(df_sp_noc[cyt], color = 'darkgreen')
	plt.legend(['Control', 'SP'])
	plt.xlabel('{} levels'.format(cyt))
	plt.savefig('./../plots/ctrl_sp/{}_noc.png'.format(cyt), dpi = 300)
	print('Saved ctrl_sp/{}'.format(cyt))
	plt.clf()

#CTRL VS RR
for cyt in cyt_list:
	plt.title('{} - RR vs CTRL(No CS Therapy)\nN = {}'.format(cyt, len(df_rr_noc)))
	sns.displot(df_ctrl_noc[cyt], color = 'grey')
	sns.displot(df_rr_noc[cyt], color = 'darkgreen')
	plt.legend(['Control', 'RR'])
	plt.xlabel('{} levels'.format(cyt))
	plt.savefig('./../plots/ctrl_rr/{}_noc.png'.format(cyt), dpi = 300)
	print('Saved ctrl_rr/{}'.format(cyt))
	plt.clf()

#creating sub-dfs for each class (therapy before PL)
df_ctrl_cort = df[(df['class_int'] == 6) & (df['cort_therapy'] == 0)].dropna() #probably empty -.-
df_pp_cort = df[(df['class_int'] == 5) & (df['cort_therapy'] == 0)].dropna()
df_sp_cort = df[(df['class_int'] == 4) & (df['cort_therapy'] == 0)].dropna()
df_rr_cort = df[(df['class_int'] == 3) & (df['cort_therapy'] == 0)].dropna()

#plotting distributions

#CTRL VS PP
for cyt in cyt_list:
	plt.title('{} - PP vs CTRL(CS Therapy)\nN = {}'.format(cyt, len(df_pp_cort)))
	sns.displot(df_ctrl_cort[cyt], color = 'grey')
	sns.displot(df_pp_cort[cyt], color = 'darkgreen')
	plt.legend(['Control', 'PP'])
	plt.xlabel('{} levels'.format(cyt))
	plt.savefig('./../plots/ctrl_pp/{}_cort.png'.format(cyt), dpi = 300)
	print('Saved ctrl_pp/{}'.format(cyt))
	plt.clf()

#CTRL VS SP
for cyt in cyt_list:
	plt.title('{} - SP vs CTRL(CS Therapy)\nN = {}'.format(cyt, len(df_sp_cort)))
	sns.displot(df_ctrl_cort[cyt], color = 'grey')
	sns.displot(df_sp_cort[cyt], color = 'darkgreen')
	plt.legend(['Control', 'SP'])
	plt.xlabel('{} levels'.format(cyt))
	plt.savefig('./../plots/ctrl_sp/{}_cort.png'.format(cyt), dpi = 300)
	print('Saved ctrl_sp/{}'.format(cyt))
	plt.clf()

#CTRL VS RR
for cyt in cyt_list:
	plt.title('{} - RR vs CTRL(CS Therapy)\nN = {}'.format(cyt, len(df_rr_cort)))
	sns.displot(df_ctrl_cort[cyt], color = 'grey')
	sns.displot(df_rr_cort[cyt], color = 'darkgreen')
	plt.legend(['Control', 'RR'])
	plt.xlabel('{} levels'.format(cyt))
	plt.savefig('./../plots/ctrl_rr/{}_cort.png'.format(cyt), dpi = 300)
	print('Saved ctrl_rr/{}'.format(cyt))
	plt.clf()