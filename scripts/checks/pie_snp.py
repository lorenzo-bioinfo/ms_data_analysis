from matplotlib import pyplot as plt
import seaborn as sns
colors = sns.color_palette('viridis')[0:2]
x = [30, 29, 44]
lab = 'G1,G2,G3'.split(',')
explode = [0.01] * 3
plt.pie(x, labels = lab, colors = colors, explode = explode)
plt.savefig('pie_snp.png', dpi = 300)
plt.clf()