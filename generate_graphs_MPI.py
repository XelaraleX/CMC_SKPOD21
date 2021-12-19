import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
import dataframe_image as dfi

if __name__ == '__main__':
    data = pd.read_csv('data.csv', sep=';')
    print(data)
    palette = sns.color_palette("flare", 6)
    sns.set()

    ax = sns.lineplot(data=data, x='n', y='time', hue='nprocs', palette=palette)
    ax.set_title('Execution time in seconds with nxnxn matrix')
    ax.set_xscale('log', basex=2)
    plt.savefig('lineplot.png')
    plt.show()

    big_data = data[data['n'] > 2 ** 6]

    ax = sns.lineplot(data=big_data, x='n', y='time', hue='nprocs', palette=palette)
    ax.set_title('Execution time in seconds with nxnxn matrix')
    ax.set_xscale('log', basex=2)
    plt.savefig('lineplot2.png')
    plt.show()

    mega_data = data[data['n'] > 2 ** 7]

    ax = sns.lineplot(data=mega_data, x='n', y='time', hue='nprocs', palette=palette)
    ax.set_title('Execution time in seconds with nxnxn matrix')
    ax.set_xscale('log', basex=2)
    plt.savefig('lineplot3.png')
    plt.show()

    data = data.pivot("n", "nprocs", "time")
    dfi.export(data, 'table.png')

    ax = sns.heatmap(data, norm=LogNorm(), square=True)
    ax.set_title('Execution time in seconds with nxnxn matrix')
    plt.savefig('heatmap.png')
    plt.show()
