import pandas as pd
from itertools import combinations
import seaborn as sns
from matplotlib import pyplot as plt

cyt_list = 'IL1B,IL2,IL4,IL5,IL6,IL7,CXCL8,IL10,IL12B,IL13,IL17A,CSF3,CSF2,IFNG,CCL2,CCL4,TNF,IL1RN,IL9,IL15,CCL11,FGF2,CXCL10,PDGFB,CCL5,VEGFA,CCL3'.split(',')
#getting all possible combinations of 2 cytokines
cyt_comb = list(combinations(cyt_list, 2))

#importing and cleaning data

df = pd.read_excel('../database/db.xlsx', sheet_name = 'SM NM', usecols = 'A,F:AF,CB,DD')
df.columns = ('patient_id,IL1B,IL2,IL4,IL5,IL6,IL7,CXCL8,IL10,IL12B,IL13,IL17A,CSF3,CSF2,IFNG,CCL2,CCL4,TNF,IL1RN,IL9,IL15,CCL11,FGF2,CXCL10,PDGFB,CCL5,VEGFA,CCL3,class_int,cort_therapy').split(',')
df_rr = df[(df['class_int'] == 3)]
df_rr_filled = df_rr.fillna({'cort_therapy' : 1})
df_rr_ok = df_rr_filled[df_rr_filled['cort_therapy'] == 1].dropna()

#plotting
for x, y in cyt_comb:
	sns.regplot(data = df_rr_ok, x = x, y = y, robust = True)
	plt.title(f'{x} v {y}')
	plt.savefig(f'../plots/cyt_regressions/{x}-{y}.png', dpi = 300)
	plt.clf()

#Note to self: by setting the option 'robust' on True, the function uses a
#more accurate version of the regression which outweights outliers.
#However.. the documentation said this was going be more computationally intensive.
#Next time, BELIEVE IT. THIS SCRIPT HAS BEEN RUNNING FOR HOURS!!