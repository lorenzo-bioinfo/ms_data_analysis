from matplotlib import pyplot as plt
import seaborn as sns
colors = sns.color_palette('mako')[0:6]
x = [36, 24, 26, 68, 38, 38]
lab = 'G1,G2,G3,G4,G5,G6'.split(',')
explode = [0.01] * 6
plt.pie(x, labels = lab, colors = colors, explode = explode)
plt.savefig('pie_norm.png', dpi = 300)
plt.clf()