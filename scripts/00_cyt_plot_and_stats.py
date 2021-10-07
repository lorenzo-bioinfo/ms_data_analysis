import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import math

#getting a list of cytokines names/labels
cyt_list = 'IL1B,IL2,IL4,IL5,IL6,IL7,CXCL8,IL10,IL12B,IL13,IL17A,CSF3,CSF2,IFNG,CCL2,CCL4,TNF,IL1RN,IL9,IL15,CCL11,FGF2,CXCL10,PDGFB,CCL5,VEGFA,CCL3'.split(',')
#getting dataframe from csv previously exported
cyt_df = pd.read_excel('../database/db.xlsx', sheet_name = 'SM NM', usecols = 'F:AF,CB')
cyt_list.append('class')
cyt_df.columns = cyt_list
cyt_list.pop()
#cleaning df from NaN values
cyt_df.dropna(inplace = True) #done: 450 rows of 490 preserved ('no liquor' out too)

#getting cyt_df for each patients' class:

cyt_ctrl = cyt_df[cyt_df['class'] == 6]
cyt_rr = cyt_df[cyt_df['class'] == 3]
cyt_pp = cyt_df[cyt_df['class'] == 5]
cyt_sp = cyt_df[cyt_df['class'] == 4]

#Getting the distribution for each cytokine and 
#superimposing it to the control cytokine distribution

sns.distributions._has_statsmodels = False #needed to avoid kde error coming from sns using statsmodel

#CTRL VS PP
for cyt in cyt_list:
	plt.title('{} - PP vs CTRL\nN = {}'.format(cyt, len(cyt_pp)))
	sns.distplot(cyt_ctrl[cyt], color = 'grey')
	sns.distplot(cyt_pp[cyt], color = 'darkgreen')
	plt.legend(['Control', 'PP'])
	plt.xlabel('{} levels'.format(cyt))
	plt.savefig('./../plots/ctrl_pp/{}.png'.format(cyt), dpi = 300)
	print('Saved ctrl_pp/{}'.format(cyt))
	plt.clf()

#CTRL VS SP
for cyt in cyt_list:
	plt.title('{} - SP vs CTRL\nN = {}'.format(cyt, len(cyt_sp)))
	sns.distplot(cyt_ctrl[cyt], color = 'grey')
	sns.distplot(cyt_sp[cyt], color = 'darkgreen')
	plt.legend(['Control', 'SP'])
	plt.xlabel('{} levels'.format(cyt))
	plt.savefig('./../plots/ctrl_sp/{}.png'.format(cyt), dpi = 300)
	print('Saved ctrl_sp/{}'.format(cyt))
	plt.clf()

#CTRL VS RR
for cyt in cyt_list:
	plt.title('{} - RR vs CTRL\nN = {}'.format(cyt, len(cyt_rr)))
	sns.distplot(cyt_ctrl[cyt], color = 'grey')
	sns.distplot(cyt_rr[cyt], color = 'darkgreen')
	plt.legend(['Control', 'RR'])
	plt.xlabel('{} levels'.format(cyt))
	plt.savefig('./../plots/ctrl_rr/{}.png'.format(cyt), dpi = 300)
	print('Saved ctrl_rr/{}'.format(cyt))
	plt.clf()

#creating dictionary for ctrl mean cytokine levels

ctrl_mean_list = []
for cyt in cyt_list:
	mean = cyt_ctrl[cyt].astype(float).mean()
	ctrl_mean_list.append(mean)

ctrl_mean_dict = dict(zip(cyt_list, ctrl_mean_list))

#getting a csv with more statistics:
cyt_lev_dfs = [cyt_ctrl, cyt_rr, cyt_pp, cyt_sp]

with open('data/cytokine_statistics/full_stats.tsv', 'w') as f:
	f.write('cytokine\tctrl_mean\tctrl_std\tpp_mean\tpp_std\tsp_mean\tsp_std\trr_mean\trr_std\tpp_diff\tsp_diff\trr_diff\nrr_d')
	for cyt in cyt_list:
		ctrl_mean = ctrl_mean_dict[cyt]
		ctrl_std = cyt_ctrl[cyt].astype(float).std()
		pp_mean = cyt_pp[cyt].astype(float).mean()
		pp_std = cyt_pp[cyt].astype(float).std()
		pp_diff = (pp_mean - ctrl_mean)/math.sqrt(pp_std * ctrl_std) #define what to do with this value
		sp_mean = cyt_sp[cyt].astype(float).mean()
		sp_std = cyt_sp[cyt].astype(float).std()
		sp_diff = (sp_mean - ctrl_mean)/math.sqrt(sp_std * ctrl_std)
		rr_mean = cyt_rr[cyt].astype(float).mean()
		rr_std = cyt_rr[cyt].astype(float).std()
		rr_diff = (rr_mean - ctrl_mean)/math.sqrt(rr_std * ctrl_std)
		rr_d = (rr_mean - ctrl_mean)
		line = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(cyt, ctrl_mean, ctrl_std, pp_mean, pp_std, sp_mean, sp_std, rr_mean, rr_std, pp_diff, sp_diff, rr_diff, rr_d)
		f.write(line)

stats_df = pd.read_csv('data/cytokine_statistics/full_stats.tsv', sep='\t')
print(stats_df)