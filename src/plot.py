import seaborn as sns
import matplotlib.pyplot as plt

def formatplot(labels, figname):
    plt.xlabel(labels['xlabel'])
    plt.ylabel(labels['ylabel'])
    plt.title(labels['title'])
    plt.savefig(figname)

def scatterplot(data, labels, figname):
    fig = plt.figure()
    sns.scatterplot(data, x=labels['xlabel'], y=labels['ylabel'])
    formatplot(labels, figname)

def heatmap(data, labels, figname):
    fig = plt.figure()
    sns.heatmap(data, x=labels['xlabel'], y=labels['ylabel'])
    formatplot(labels, figname)

def barchart(data, labels, figname):
    fig = plt.figure()
    sns.barchart(data, x=labels['xlabel'], y=labels['ylabel'])
    formatplot(labels, figname)
