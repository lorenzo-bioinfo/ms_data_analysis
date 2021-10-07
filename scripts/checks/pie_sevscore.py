from matplotlib import pyplot as plt
import seaborn as sns
colors = sns.color_palette('summer')[0:3]
x = [38, 61, 51, 60]
lab = 'G1,G2,G3,G4'.split(',')
explode = [0.01] * 4
plt.pie(x, labels = lab, colors = colors, explode = explode)
plt.savefig('pie_sevscores.png', dpi = 300)
plt.clf()