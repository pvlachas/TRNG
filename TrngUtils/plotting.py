#!/usr/bin/env python
# # -*- coding: utf-8 -*-

"""Created by: Vlachas Pantelis, CSE-lab, ETH Zurich
"""
#!/usr/bin/env python


import numpy as np
from statsmodels.graphics.tsaplots import plot_acf

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)
rc('font', family='serif')

def makeAutocorrelationPlot(dataseries, lags, dt, T_l, path):
    assert len(np.shape(dataseries))==1

    fig, ax = plt.subplots()
    plot_acf(dataseries, lags=lags, ax=ax, title=None)
    line = ax.lines[1]
    autocorrelation = line.get_ydata()

    tmax = lags * dt
    nmax = int(np.floor(tmax/T_l))
    newxaxis = np.arange(nmax) * T_l / dt
    # Adding additional axis with Lyapunov time
    new_tick_locations = np.array(newxaxis)
    def tick_function(X):
        V = X * dt / T_l
        return ["%.3f" % z for z in V]
    ax.set_xlim([0, lags])
    ax.set_xlabel("Lag")
    ax.set_ylabel("Autocorrelation")
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(new_tick_locations)
    ax2.set_xticklabels(tick_function(new_tick_locations))
    ax2.set_xlabel(r"$t \, / \, \Lambda_1^{-1}$")
    # plt.show()
    fig.savefig(path)
    plt.close()

    return autocorrelation


def makeHistogramPlot(data, nbins, figure_path):
    plt.hist(data, bins = nbins)
    plt.xlabel("Value")
    plt.ylabel("Counts")
    plt.savefig(figure_path)
    plt.close()

def plotECDF(x, y, ecdf_func, figure_path):
    plt.plot(x, y, drawstyle='steps-post', label='ecdf')
    plt.plot(x, ecdf_func(x), label='fit')
    plt.xlabel("$x$")
    plt.ylabel("$F_X(x)$")
    plt.grid(True)
    plt.savefig(figure_path)
    plt.close()










