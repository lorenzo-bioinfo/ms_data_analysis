import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import os


df = pd.read_excel('../database/db.xlsx', sheet_name = 'SM NM', usecols = 'A,F:AF,CB,CK,CM')
df.columns = ('patient_id,IL1B,IL2,IL4,IL5,IL6,IL7,CXCL8,IL10,IL12B,IL13,IL17A,CSF3,CSF2,IFNG,CCL2,CCL4,TNF,IL1RN,IL9,IL15,CCL11,FGF2,CXCL10,PDGFB,CCL5,VEGFA,CCL3,class_int,gender,age').split(',')
df.dropna(inplace = True)
df.to_csv('data/cyt_age_df.csv')

#gender column values: 1 = F, 0 = m

#searching minimum and maximum age in df

younger = min(df['age'].astype(float))
older = max(df['age'].astype(float))
print('Min age is: ', younger)
print('Max age is: ', older)

#defining ranges of age in a list

age_range = [(float(younger), 27), (27, 37), (37, 47), (47, 57), (57, 67), (67, float(older))]

#creating sub-dfs for age ranges

def getAgeDf(df, interval):
	age_df = df[(df['age'] >= interval[0]) & (df['age'] < interval[1])]
	return age_df

df_1627 = getAgeDf(df, age_range[0])
df_2737 = getAgeDf(df, age_range[1])
df_3747 = getAgeDf(df, age_range[2])
df_4757 = getAgeDf(df, age_range[3])
df_5767 = getAgeDf(df, age_range[4])
df_6777 = getAgeDf(df, age_range[5])

age_dfs = [df_1627, df_2737, df_3747, df_4757, df_5767, df_6777]
legend_labels = ['16-27', '27-37', '37-47', '45-47', '57-67', '67-77']
legend_colors = 'blue,darkgreen,red,black,purple,navy'.split(',')
#all categories are quite similar in number, except last two
	
#plotting distributions for each age range:	

cyt_list = 'IL1B,IL2,IL4,IL5,IL6,IL7,CXCL8,IL10,IL12B,IL13,IL17A,CSF3,CSF2,IFNG,CCL2,CCL4,TNF,IL1RN,IL9,IL15,CCL11,FGF2,CXCL10,PDGFB,CCL5,VEGFA,CCL3'.split(',')
sns.distributions._has_statsmodels = False #needed to avoid kde error coming from sns using statsmodel

for cyt in cyt_list:
	i = 0
	for dataf in age_dfs:
		sns.kdeplot(dataf[cyt].astype(float), color = legend_colors[i])
		i += 1
	plt.title('{} - Age comparison'.format(cyt))
	plt.legend(legend_labels)
	filename = '../plots/cyt_age/{}'.format(cyt)
	plt.savefig(filename, dpi = 300)
	print('Saved plot {}'.format(filename))
	plt.clf()