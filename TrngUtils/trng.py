#!/usr/bin/env python
# # -*- coding: utf-8 -*-

"""Created by: Vlachas Pantelis, CSE-lab, ETH Zurich
"""
#!/usr/bin/env python

import math
import numpy as np
from scipy.interpolate import interp1d

from . import plotting
from . import utils

def ecdf(samples):
    # Computes empirical CDF from samples
    x, counts = np.unique(samples, return_counts=True)
    cusum = np.cumsum(counts)
    return x, cusum / cusum[-1]

def computeECDF(data):
    assert len(np.shape(data))==1
    x, y = ecdf(data)
    ecdf_func = interp1d(x, y, kind='cubic', fill_value="extrapolate")
    return ecdf_func, x, y

def performDecorrelationAnalysis(u, dt, lags, T_l, saving_path_figures, saving_path_data):
    # Total number of dimensions
    D = np.shape(u)[1]
    for proj in range(D):
        dataseries = u[:, proj]
        fig_path = saving_path_figures + '/Autocorrelation_plot_proj{:}.pdf'.format(proj)
        autocorrelation = plotting.makeAutocorrelationPlot(dataseries, lags, dt, T_l, fig_path)
        zero_crossings = np.where(np.diff(np.sign(autocorrelation)))[0]
        print(zero_crossings)
        # Decorrelation time
        T_d = zero_crossings[0] * dt
        print("Lyapunov time= {:}".format(T_l))
        print("Decorrelation time= {:}".format(T_d))

        data_dict = {
        "u":u,
        "proj":proj,
        "T_l":T_l,
        "T_d":T_d,
        "zero_crossings":zero_crossings,
        "autocorrelation":autocorrelation,
        }
        data_path = saving_path_data + '/Autocorrelation_proj{:}.pickle'.format(proj)
        utils.saveData(data_dict, data_path)
    return 0


def runsTest(data, data_median):
    runs, n1, n2 = 0, 0, 0
    # Checking for start of new run
    for i in range(len(data)):
          
        # no. of runs
        if (data[i] >= data_median and data[i-1] < data_median) or \
                (data[i] < data_median and data[i-1] >= data_median):
            runs += 1  
          
        # no. of positive values
        if(data[i]) >= data_median:
            n1 += 1   
          
        # no. of negative values
        else:
            n2 += 1   
  
    runs_exp = ((2*n1*n2)/(n1+n2))+1
    stan_dev = math.sqrt((2*n1*n2*(2*n1*n2-n1-n2))/ \
                       (((n1+n2)**2)*(n1+n2-1)))
  
    z = (runs-runs_exp)/stan_dev
  
    return z







