import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

if __name__ == '__main__':
    data = pd.read_csv('data.csv', sep='\t')

    palette = sns.color_palette("flare", 7)
    sns.set()

    ax = sns.lineplot(data=data, x='n', y='time', hue='threads', palette=palette)
    ax.set_title('Execution time in seconds with nxnxn matrix')
    ax.set_xscale('log', basex=2)
    plt.savefig('lineplot.png')
    plt.show()

    data = data.pivot("n", "threads", "time")

    ax = sns.heatmap(data, norm=LogNorm(), square=True)
    ax.set_title('Execution time in seconds with nxnxn matrix')
    plt.savefig('heatmap.png')
    plt.show()
