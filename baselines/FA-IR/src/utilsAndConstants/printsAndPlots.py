'''
Created on Apr 12, 2017

@author: meike.zehlike
'''
import numpy as np
import matplotlib as mpl

from matplotlib import pyplot as plt


def plotFourListsInOnePlot(xdata, ydata1, ydata2, ydata3, ydata4, xlabel, ylabel, filename):
    mpl.rcParams.update({'font.size': 20, 'lines.linewidth': 3, 'lines.markersize': 15})
    # avoid type 3 (i.e. bitmap) fonts in figures
    mpl.rcParams['ps.useafm'] = True
    mpl.rcParams['pdf.use14corefonts'] = True
    mpl.rcParams['text.usetex'] = True
    plt.plot(xdata, ydata1, c='r')
    plt.scatter(xdata, ydata2, marker='x', c='r', s=100)
    plt.plot(xdata, ydata3, c='b')
    plt.scatter(xdata, ydata4, marker='x', c='b', s=100)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid()
    plt.tight_layout()

    if filename:
        plt.savefig(filename)
    else:
        plt.show()

def plotTwoListsInOnePlot(xdata, ydata1, ydata2, labelLine1, labelLine2,
                          labelXAx,
                          labelYAx1,
                          labelYAx2,
                          filename):
    mpl.rcParams.update({'font.size': 20, 'lines.linewidth': 3, 'lines.markersize': 15, 'font.family':'Times New Roman'})
    # avoid type 3 (i.e. bitmap) fonts in figures
    mpl.rcParams['ps.useafm'] = True
    mpl.rcParams['pdf.use14corefonts'] = True
    mpl.rcParams['text.usetex'] = True
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax2 = ax.twinx()
    ax.plot(xdata, ydata1, marker='.', c='r', label=labelLine1)
    ax2.plot(xdata, ydata2, linestyle='--', marker='.', c='b', label=labelLine2)

    # added these two lines
    ax.legend(loc=1, framealpha=0.7)
    ax2.legend(loc=4, framealpha=0.7)

    ax.set_xlim(xmin=0, xmax=1)
    ax.set_ylim(ymin=0, ymax=1)

    ax.set_xlabel(labelXAx)
    ax.set_ylabel(labelYAx1)
    ax2.set_ylabel(labelYAx2)

    fig.tight_layout()
    plt.savefig(filename)
