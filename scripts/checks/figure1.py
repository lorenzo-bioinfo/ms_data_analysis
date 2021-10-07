import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
cyt_list = 'IL1B,IL2,IL4,IL5,IL6,IL7,CXCL8,IL10,IL12B,IL13,IL17A,CSF3,CSF2,IFNG,CCL2,CCL4,TNF,IL1RN,IL9,IL15,CCL11,FGF2,CXCL10,PDGFB,CCL5,VEGFA,CCL3'.split(',')
df = pd.read_excel('../../database/db.xlsx', sheet_name = 'SM NM', usecols = 'A,F:AF,CB,CN,DD')
df.columns = 'patient_id,IL1B,IL2,IL4,IL5,IL6,IL7,CXCL8,IL10,IL12B,IL13,IL17A,CSF3,CSF2,IFNG,CCL2,CCL4,TNF,IL1RN,IL9,IL15,CCL11,FGF2,CXCL10,PDGFB,CCL5,VEGFA,CCL3,class_int,act,cort_therapy'.split(',')
df_rr = df[df['class_int'] == 3]
df_rr_filled = df_rr.fillna({'cort_therapy' : 1, 'act' : 0})
df_ok = df_rr_filled[df_rr_filled['cort_therapy'] == 1].dropna()

for cyt in cyt_list:
	vect = list(df_ok[cyt].astype(float))
	sns.histplot(vect, stat = 'density', kde = True, color = 'darkgreen')
	plt.title(f'{cyt}')
	plt.xlabel('CSF Levels')
	plt.ylabel('Frequency')
	plt.savefig(f'./fig1/{cyt}.png', dpi = 300)
	plt.clf()