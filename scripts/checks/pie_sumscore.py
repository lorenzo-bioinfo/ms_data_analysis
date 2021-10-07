from matplotlib import pyplot as plt
import seaborn as sns
colors = sns.color_palette('summer')[0:4]
x = [27, 49, 50, 58, 46]
lab = 'G1,G2,G3,G4,G5'.split(',')
explode = [0.01] * 5
plt.pie(x, labels = lab, colors = colors, explode = explode)
plt.savefig('pie_sumscore.png', dpi = 300)
plt.clf()