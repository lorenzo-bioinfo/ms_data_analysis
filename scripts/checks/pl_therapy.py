import pandas as pd
import datetime as dt

with open('../data/rr_noc_list.txt', 'r') as f:
	pat_list = f.readline().strip().split(',')

df = pd.read_excel('../../database/db.xlsx', sheet_name = 'SM NM', usecols = 'A,CA,DX,DY')
cols = ['patient_id', 'pl_date', 'ther_start', 'ther1']
df.columns = cols

df_rr = df[df['patient_id'].astype(str).isin(pat_list)].dropna()
print(df_rr.columns)


df_rr.drop(df_rr.loc[df_rr['ther_start'] == 'NO'].index, inplace = True)

fl_ther = [1, 2, 3, 4, 5, 6, 12,13, 14, 15, 16, 17, 18, 19, 21, 22, 24, 24, 31]
terapie = list(df_rr['ther1'].astype(int))
linea_ther = ['FL' if x in fl_ther else 'SL' for x in terapie]
df_rr['TL'] = linea_ther
print(df_rr)
rr_list = list(df_rr['patient_id'])
pldate = list(df_rr['pl_date'])
therstart = list(df_rr['ther_start'])
days = []
for i in range(len(pldate)):
	delta = therstart[i] - pldate[i]
	days.append(delta.days)
usable_days = [x for x in days if abs(x) < 90]

df_rr['delta_pl'] = days
df_ok = df_rr[abs(df_rr['delta_pl']) < 91]

# Adding therapy 1 end date

df2 = pd.read_excel('../../database/db.xlsx', sheet_name = 'SM NM', usecols = 'A,DZ,EA,EB,EC')
columns2 = 'patient_id,end1,reason1,start2,ther2'.split(',')
df2.columns = columns2
df_rr2 = df2[df2['patient_id'].isin(rr_list)]
end1 = list(df_rr2['end1'].astype(str))
print(end1)
end1_clean = [None if x == 'NaT' else x for x in end1]
print(end1_clean)
df_rr['end1'] = end1_clean
print(df_rr)

# Adding therapy 1 end reason

r1 = list(df_rr2['reason1'].astype(str))
r1_clean = []
for el in r1:
	if el == '1.0':
		r1_clean.append('A')
	elif el == '0.0':
		r1_clean.append('P')
	else:
		r1_clean.append(None)
print(r1_clean, len(r1_clean))
df_rr['re_end1'] = r1_clean

#adding therapy 2 start date and difference in days between
#starting day of therapy 2 and ending day of therapy 1

ed1 = list(df_rr2['end1'])
sd2 = list(df_rr2['start2'])
print(ed1, sd2)

ed1_clean = []
sd2_clean = []

for x in ed1:
	if str(x) == 'NaT':
		ed1_clean.append(None)
	else:
		ed1_clean.append(str(x).split(' ')[0])
print(ed1_clean)


for x in sd2:
	if str(x) == 'NaT':
		sd2_clean.append(None)
	else:
		sd2_clean.append(str(x).split(' ')[0])
print(sd2_clean)

delta1_2 = []

for i in range(len(ed1)):
	try:
		a = sd2[i] - ed1[i]
		delta1_2.append(a.days)
	except ValueError:
		delta1_2.append(None)

print(delta1_2)
df_rr['end1'] = ed1_clean
df_rr['delta1_2'] = delta1_2
df_rr['ther2_start'] = sd2_clean

'''
for i in range(len(ed1_clean)):
	try:
		a = sd2_clean[i] - ed1_clean[i]
		delta1_2.append(a.days)
		print('No TypeError', sd2_clean[i], ed1_clean[i])
	except TypeError:
		delta1_2.append(None)
'''

df_rr.to_csv('delta_pl_ther.tsv', sep = '\t')