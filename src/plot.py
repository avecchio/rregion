import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import sys
sys.setrecursionlimit(50000)

def formatplot(labels, figname):
    plt.xlabel(labels['xlabel'])
    plt.ylabel(labels['ylabel'])
    plt.title(labels['title'])
    plt.savefig(figname)
    plt.close()

def boxen(data, figname, labels):
    fig = plt.figure()
    ax = sns.boxenplot(data, x=labels['xdata'], y=labels['hdata'], gap=.2)
    #if (labels['xlogscale']):
    #    plt.xscale('log')
    formatplot(labels, figname)

def lineplot(data, figname, labels):
    fig = plt.figure(figsize=(12, 6))
    ax = sns.lineplot(data, x=labels['xdata'], y=labels['ydata'], hue=labels['hue'])
    #plt.yscale('log')
    #plt.tight_layout()
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1,1))
    formatplot(labels, figname)

def displot(data, figname, labels):
    fig = plt.figure()
    sns.displot(data, x=labels['xdata'], hue=labels['hdata'], kind="kde")
    if (labels['xlogscale']):
        plt.xscale('log')
    if (labels['ylogscale']):
        plt.yscale('log')
    formatplot(labels, figname)

def violin(data, figname, labels):
    fig = plt.figure()
    sns.violinplot(data, x=labels['xdata'], y=labels['hdata'])
    plt.xscale('log')
    formatplot(labels, figname)

def histogram(data, figname, labels):
    fig = plt.figure()
    sns.histplot(data, x=labels['xdata'], hue=labels['hdata'], multiple="stack")
    formatplot(labels, figname)

def scatterplot(data, figname, labels, datalabels=None):
    #print(data.shape())
    plt.clf()
    print(figname, data.size, data.shape, data.columns)
    fig = plt.figure()
    if (labels['ylogscale'] or labels['xlogscale']) and (datalabels is not None):
        sns.scatterplot(data, x=labels['xdata'], y=labels['ydata'], hue=datalabels, alpha=0.3)
    elif datalabels is not None:
        sns.scatterplot(data, x=labels['xdata'], y=labels['ydata'], hue=datalabels, alpha=0.3)
    else:
        sns.scatterplot(data, x=labels['xdata'], y=labels['ydata'])
    if (labels['xlogscale']):
        plt.xscale('log')
    if (labels['ylogscale']):
        plt.yscale('log')
    formatplot(labels, figname)

def clustermap(data, figname, labels):
    fig = plt.figure(figsize=(8, 6))
    g = sns.clustermap(data)
    ax = g.ax_heatmap
    #print(ax)
    ax.set_xlabel(labels['xlabel'])
    ax.set_ylabel(labels['ylabel'])
    plt.suptitle(labels['title']) #, y=0.95)
    plt.tight_layout()
    plt.savefig(figname)
    plt.close()

def heatmap(data, figname, labels):
    fig = plt.figure()
    if 'categorize' in labels:
        categorical_column = labels['categorize']
        num_unique_categories = len(data[categorical_column].unique())
        #print('unique categories', num_unique_categories)
        color_pal = sns.color_palette('pastel', n_colors=num_unique_categories)
        #print(len(color_pal))
        colors_lut = dict(zip(list(data[categorical_column].unique()), color_pal))
        #print(colors_lut)

        categorical_colors = data[categorical_column].map(colors_lut)
        data = data.drop(categorical_column, axis=1)
        #print(data.columns)
        #print(list(categorical_colors)[0:3])
        #cbar_pos = (0.95, 0.75, 0.03, 0.15)
        g = sns.clustermap(data, row_colors=categorical_colors, yticklabels=False) #, cbar_pos=cbar_pos) #, row_cluster=False, col_cluster=False)
        #plt.legend(loc='upper right')
        #g.set(xlabel=labels['xlabel'], ylabel=labels['ylabel'])
        ax = g.ax_heatmap
        #print(ax)
        ax.set_xlabel(labels['xlabel'])
        ax.set_ylabel(labels['ylabel'])

        legend_elements = [plt.Line2D([0], [0], marker=0, color=color, markersize=10, label=label, markerfacecolor=color) for label, color in colors_lut.items()]
        legend = plt.legend(handles=legend_elements, title='Region Types', loc='upper left')
        plt.gca().add_artist(legend)
        #cbar = plt.colorbar()
        #plt.xlabel(labels['xlabel'])
        #plt.ylabel(labels['ylabel'])
        #ax.set_title(labels['title'], pad=50, loc='center')
        #for label in num_unique_categories:
        #    g.ax_col_dendrogram.bar
        plt.suptitle(labels['title']) #, y=0.95)
        plt.tight_layout()
        plt.savefig(figname)
        plt.close()
    else:
        sns.heatmap(data)
        #plt.tight_layout()
        formatplot(labels, figname)

def barchart(data, labels, figname):
    fig = plt.figure()
    sns.barchart(data, x=labels['xdata'], y=labels['ydata'])
    formatplot(labels, figname)
