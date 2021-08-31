#!/usr/bin/env python
# # -*- coding: utf-8 -*-

"""Created by: Vlachas Pantelis, CSE-lab, ETH Zurich
"""
#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import argparse


import numpy as np
import random
import statistics

from scipy.interpolate import interp1d

sys.path.append('..')
from TrngUtils import trng
from TrngUtils import plotting
from TrngUtils import utils

from Utils import chua_circut


data_path   = "../../Data/ChuaSim/Data/testing_data_N100000.pickle"
data_field  = "test_input_sequence"

saving_path_figures = "./Figures"
os.makedirs(saving_path_figures, exist_ok=True)

saving_path_data = "./Data"
os.makedirs(saving_path_data, exist_ok=True)

data = utils.loadData(data_path)
u = data[data_field]
dt = data["dt"]


print("Shape of data = {:}".format(np.shape(u)))
print("Timestep dt = {:}".format(dt))
min_u = np.min(u)
max_u = np.max(u)

N_plot = 100000
# Total number of dimensions
D = np.shape(u)[1]

for proj in range(D):
    dataseries = u[:N_plot, proj]
    lags = 100
    # Lyapunov time
    T_l = 4.17
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



print(ark)

# Building the true random number generator
trng_lag = 20

trng_dataseries = dataseries[::trng_lag]

print(np.shape(trng_dataseries))

plt.hist(trng_dataseries, bins = 20)
plt.show()


x, y = ecdf(trng_dataseries)
# x = np.insert(x, 0, x[0])
# y = np.insert(y, 0, 0.)
# print(np.any(x[1:] <= x[:-1]))

# print(ark)
ecdf_func = interp1d(x, y, kind='cubic')

plt.plot(x, y, drawstyle='steps-post', label='ecdf')
plt.plot(x, ecdf_func(x), label='fit')
plt.grid(True)
# plt.savefig('ecdf.png')
plt.show()


# Sample a new number
print(np.min(x))
print(np.max(x))
xnew = 0.4524132
sample = ecdf_func(xnew)
print(sample)

print(np.shape(trng_dataseries))
samples = np.array([ecdf_func(xnew) for xnew in trng_dataseries])
print(samples)
plt.hist(samples, bins = 100)
plt.show()


    
# Making a list of 100 random numbers 
l = []
for i in range(100):
    l.append(random.random())
      
samples_median= np.median(samples)
  
Z = abs(runsTest(samples, samples_median))
  
print('Z-statistic= ', Z)
Zcritical = 1.96 # for confidence level of 95%
print("Confidence interval 95%, $Z_{crit}="+"{:}".format(Zcritical)+"$")
if Z > Zcritical:
    print("The null hypothesis is rejected, i.e. the numbers are declared not to be random")
else:
    print("The null hypothesis is accepted, i.e. the numbers are declared random.")

# Compare the value of the calculated Z-statistic with Zcritical  for a given level of confidence (Zcritical =1.96 for confidence level of 95%) . The null hypothesis is rejected i.e. the numbers are declared not to be random, if |Z|>Zcritical . 


print(ark)










