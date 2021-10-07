from matplotlib import pyplot as plt
import seaborn as sns
colors = sns.color_palette('viridis')[0:6]
x = [36, 36, 40, 48, 24, 16, 30]
lab = 'G1,G2,G3,G4,G5,G6,G7'.split(',')
explode = [0.01] * 7
plt.pie(x, labels = lab, colors = colors, explode = explode)
plt.savefig('pie_quart.png', dpi = 300)
plt.clf()