#trying to see (graphically) if there is any degree of covariance between
#cytokines in the same group. The groups are obtained by clustering based
#on interactors

import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

cyt_groups = [['CXCL8', 'CCL4', 'CXCL10', 'CCL5'], ['IL12B', 'PDGFB', 'FGF2', 'VEGFA'], ['IFNG', 'IL2', 'IL7', 'IL15', 'CSF3', 'CSF2', 'IL10', 'IL17A', 'IL4', 'IL9', 'IL5', 'IL13'], ['TNF', 'IL6', 'IL1B', 'IL1RN', 'CCL3', 'CCL2', 'CCL11']]

#reading and cleaning data
#this time I will be using all the data, both from controls and patients across
#all classes, only excluding the ones who underwent cortisonic therapy before PL
df = pd.read_excel('../database/db.xlsx', sheet_name = 'SM NM', usecols = 'A,F:AF,CB,DD')
df.columns = ('patient_id,IL1B,IL2,IL4,IL5,IL6,IL7,CXCL8,IL10,IL12B,IL13,IL17A,CSF3,CSF2,IFNG,CCL2,CCL4,TNF,IL1RN,IL9,IL15,CCL11,FGF2,CXCL10,PDGFB,CCL5,VEGFA,CCL3,class_int,cort_therapy').split(',')
df_filled = df.fillna({'cort_therapy' : 1})
df_ok = df_filled[df_filled['cort_therapy'] == 1].dropna()
df_ok['class_int'].replace(to_replace = [1, 2, 3, 4, 5, 6, 7, 8, 9], value = ['RIS', 'CIS', 'RR', 'SP', 'PP', 'CTRL_NOINF', 'CTRL_INF', 'ALTRO', 'GILENYA'], inplace = True)
classes = ['RR', 'CTRL_NOINF']
df_clean = df_ok[df_ok['class_int'].isin(classes)]
#plotting

for idx, group in enumerate(cyt_groups):
	sns.pairplot(df_clean, hue = 'class_int', vars = group, kind = 'reg', diag_kind = 'hist')
	plt.savefig(f'./../plots/cyt_pairplots/group{idx + 1}.png', dpi = 300)
	plt.title(f'Group {idx + 1} cytokines')
	plt.clf()