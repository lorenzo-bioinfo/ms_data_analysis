import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
cyt_list = 'IL1B,IL2,IL4,IL5,IL6,IL7,CXCL8,IL10,IL12B,IL13,IL17A,CSF3,CSF2,IFNG,CCL2,CCL4,TNF,IL1RN,IL9,IL15,CCL11,FGF2,CXCL10,PDGFB,CCL5,VEGFA,CCL3'.split(',')
df = pd.read_excel('../../database/db.xlsx', sheet_name = 'SM NM', usecols = 'A,F:AF,CB,CN,DD')
df.columns = 'patient_id,IL1B,IL2,IL4,IL5,IL6,IL7,CXCL8,IL10,IL12B,IL13,IL17A,CSF3,CSF2,IFNG,CCL2,CCL4,TNF,IL1RN,IL9,IL15,CCL11,FGF2,CXCL10,PDGFB,CCL5,VEGFA,CCL3,class_int,act,cort_therapy'.split(',')
df_rr = df[df['class_int'] == 3]
df_rr_filled = df_rr.fillna({'cort_therapy' : 1, 'act' : 0})
df_ok = df_rr_filled[df_rr_filled['cort_therapy'] == 1].dropna()
df_ct = df[df['class_int'] == 6]
df_ct_filled = df_ct.fillna({'cort_therapy' : 1, 'act' : 0})
df_ok_ct = df_ct_filled[df_ct_filled['cort_therapy'] == 1].dropna()

for cyt in cyt_list:
	vect1 = list(df_ok[cyt].astype(float))
	vect2 = list(df_ok_ct[cyt].astype(float))
	sns.histplot(vect1,stat = 'density', kde = True, color = 'purple')
	sns.histplot(vect2,stat = 'density', kde = True, color = 'darkgreen' )
	plt.title(f'{cyt}')
	plt.xlabel('CSF Levels')
	plt.ylabel('Frequency')
	plt.legend(['RR', 'CTRLS'])
	plt.savefig(f'./fig4/{cyt}.png', dpi = 300)
	plt.clf()

for cyt in cyt_list:
	ser1 = df_ok[cyt].astype(float)
	qnt1 = ser1.quantile(.95)
	vect1 = list(ser1.drop(ser1.index[ser1 > qnt1]))
	ser2 = df_ok_ct[cyt].astype(float)
	qnt2 = ser2.quantile(.95)
	vect2 = list(ser2.drop(ser2.index[ser2 > qnt2]))
	sns.histplot(vect1,stat = 'density', kde = True, color = 'purple')
	sns.histplot(vect2,stat = 'density', kde = True, color = 'darkgreen' )
	plt.title(f'{cyt} - No Outliers')
	plt.xlabel('CSF Levels')
	plt.ylabel('Frequency')
	plt.legend(['RR', 'CTRLS'])
	plt.savefig(f'./fig5/{cyt}.png', dpi = 300)
	plt.clf()
