import seaborn as sns
import matplotlib.pyplot as plt

def formatplot(labels, figname):
    plt.xlabel(labels['xlabel'])
    plt.ylabel(labels['ylabel'])
    plt.title(labels['title'])
    plt.savefig(figname)

def boxen(data, figname, labels):
    fig = plt.figure()
    sns.boxenplot(data, x=labels['xdata'], y=labels['hdata'], gap=.2)
    if (labels['xlogscale']):
        plt.xscale('log')
    #if (labels['ylogscale']):
    #    plt.yscale('log')
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

def scatterplot(data, figname, labels):
    fig = plt.figure()
    sns.scatterplot(data, x=labels['xdata'], y=labels['ydata'])
    formatplot(labels, figname)

def heatmap(data, labels, figname):
    fig = plt.figure()
    sns.heatmap(data, x=labels['xdata'], y=labels['ydata'])
    formatplot(labels, figname)

def barchart(data, labels, figname):
    fig = plt.figure()
    sns.barchart(data, x=labels['xdata'], y=labels['ydata'])
    formatplot(labels, figname)
