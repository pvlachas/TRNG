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


sys.path.append('..')
from TrngUtils import trng
from TrngUtils import plotting
from TrngUtils import utils

from sp800_22_tests import *

from Utils import chua_circut

Ν = 1879019
data_path   = "../../Data/ChuaSim/Data/testing_data_N{:}.pickle".format(Ν)
data_field  = "test_input_sequence"

saving_path_figures = "./Figures"
os.makedirs(saving_path_figures, exist_ok=True)

saving_path_data = "./Data"
os.makedirs(saving_path_data, exist_ok=True)

data = utils.loadData(data_path)
dataseries = data[data_field]
dt = data["dt"]


print("Shape of data = {:}".format(np.shape(dataseries)))
print("Timestep dt = {:}".format(dt))

N_plot = 1000000000
dataseries = dataseries[:N_plot]
D = np.shape(dataseries)[1]

# lags = 100
# # Lyapunov time
# T_l = 4.17
# trng.performDecorrelationAnalysis(dataseries, dt, lags, T_l, saving_path_figures, saving_path_data)




"""
Building a TRNG based on
1. A projection (x, y, or z) 
2. A time-lag
"""

# proj = 0
# trng_lag = 20



# lagtimes = [1] + list(range(0, 1000, 50))[1:]
lagtimes = list(range(100, 2000, 50))


# for proj in range(D):
for proj in [0]:
    for trng_lag in lagtimes:
    # for trng_lag in [850]:
        print("-"*80)
        print("Projection {:}, lag-time {:}".format(proj, trng_lag))
        trng_data = dataseries[::trng_lag]
        trng_dataseries = trng_data[:, proj]
        print("Shape of dataseries used in TRNG: {:}".format(np.shape(trng_dataseries)))

        nsamples = len(trng_dataseries)
        nbins = utils.getNumberOfBins(nsamples, rule="rice")

        figure_path = saving_path_figures + "/Histogram_proj{:}_lag{:}.pdf".format(proj, trng_lag)
        # plotting.makeHistogramPlot(trng_dataseries, nbins, figure_path)

        """ Estimate the empirical cummulative distribution function """
        ecdf_func, ecdf_x, ecdf_y = trng.computeECDF(trng_dataseries)
        figure_path = saving_path_figures + "/ECDF_proj{:}_lag{:}.pdf".format(proj, trng_lag)
        # plotting.plotECDF(ecdf_x, ecdf_y, ecdf_func, figure_path)

        # """
        # Example on how to get a sample form uniform,
        # given a sample from the chaotic time-series
        # """
        # # Sample a datapoint from the chaotic time-series
        # print(np.min(ecdf_x))
        # print(np.max(ecdf_x))
        # xnew = 0.4524132
        # sample = ecdf_func(xnew)
        # print(sample)


        """ Check that the samples are uniform """
        samples_uniform = np.array([ecdf_func(xnew) for xnew in trng_dataseries])
        figure_path = saving_path_figures + "/Histogram_proj{:}_lag{:}_uniform.pdf".format(proj, trng_lag)
        nbins = utils.getNumberOfBins(len(samples_uniform), rule="rice")
        # plotting.makeHistogramPlot(samples_uniform, nbins, figure_path)

        samples_median = np.median(samples_uniform)

        """ Write samples to file """
        filename_binseq = saving_path_data + "/BinaryRandSequence_proj{:}_lag{:}.bin".format(proj, trng_lag)
        # print(samples_uniform)
        # print(samples_median)

        samples_binary = np.zeros_like(samples_uniform)

        samples_binary[samples_uniform>=samples_median] = 1
        samples_binary[samples_uniform<samples_median] = 0

        utils.writeDataToByteArray(samples_binary, filename_binseq)

        # print(dir(sp800_22_tests))
        num_tests_passed, total_tests, success_rate, results = sp800_22.runTestRoutines(filename_binseq)

        # print(ark)
        # Z = np.abs(trng.runsTest(samples_uniform, samples_median))
        # print('Z-statistic= ', Z)
        # Zcritical = 1.96 # for confidence level of 95%
        # print("Confidence interval 95%, $Z_{crit}="+"{:}".format(Zcritical)+"$")
        # if Z > Zcritical:
        #     print("The null hypothesis is rejected, i.e. the numbers are declared not to be random")
        # else:
        #     print("The null hypothesis is accepted, i.e. the numbers are declared random.")

        # """
        # Compare the value of the calculated Z-statistic with Zcritical
        # for a given level of confidence (Zcritical =1.96 for confidence level of 95%).
        # The null hypothesis is rejected i.e. the numbers are declared not to be random, if |Z|>Zcritical.
        # """

        data = {
            "num_tests_passed":num_tests_passed,
            "total_tests":total_tests,
            "success_rate":success_rate,
            "results":results,
            "trng_data":trng_data,
            "trng_dataseries":trng_dataseries,
            "ecdf_func":ecdf_func,
            "ecdf_x":ecdf_x,
            "ecdf_y":ecdf_y,
            "samples_uniform":samples_uniform,
            "samples_median":samples_median,
            # "Z":Z,
        }
        data_path = saving_path_data + '/RandomNumberTests_proj{:}_lag{:}.pickle'.format(proj, trng_lag)
        utils.saveData(data, data_path, add_file_format=False)





# import matplotlib as mpl
# import matplotlib.pyplot as plt
# from matplotlib import rc
# rc('text', usetex=True)
# rc('font', family='serif')

# Zcritical = 1.96

# for proj in range(D):
#     z_values = []
#     # for trng_lag in [1, 20, 100]:
#     # for trng_lag in [250]:
#     for trng_lag in lagtimes:

#         data_path = saving_path_data + '/RandomNumberTests_proj{:}_lag{:}.pickle'.format(proj, trng_lag)
#         data = utils.loadData(data_path, add_file_format=False)
#         Z = data["Z"]
#         z_values.append(Z)

#     plt.plot(lagtimes, z_values)

# plt.plot(lagtimes, Zcritical * np.ones_like(lagtimes), "r--")
# plt.ylim([0, 2*Zcritical])
# plt.show()




"""
How much noise can we afford ?
If the system state can be measured very accurately, since the dynamics are known, the random numbers can be predicted with very high correlation.
Adding noise to the system deteriorates the situation.
"""




# proj = 1
# for trng_lag in lagtimes:

#     data_path = saving_path_data + '/RandomNumberTests_proj{:}_lag{:}.pickle'.format(proj, trng_lag)
#     data = utils.loadData(data_path, add_file_format=False)

#     trng_data = data["trng_data"]
#     samples_uniform = data["samples_uniform"]
#     ecdf_func = data["ecdf_func"]

#     """ Propagate the dynamics forward """
#     t_f = trng_lag * dt

#     for noise_level in [0, 1, 5, 10]:
#         print("-"*80)
#         print("Projection {:}, lag-time {:}, noise level {:}".format(proj, trng_lag, noise_level))
#         """ Add noise to the system """
#         data_std = np.std(trng_data, axis=0)
#         noise_var = np.power(noise_level/100. * data_std, 2.0)
#         # print(np.var(trng_data, axis=0))
#         # print(noise_var)
#         trng_data_noisy = [u0 + noise_var * np.random.randn() for u0 in trng_data]
#         trng_data_recon = np.array([chua_circut.propagateChuaDynamics(u0, t_f) for u0 in trng_data_noisy])
#         # print(np.shape(trng_data_recon))
#         trng_dataseries_recon = trng_data_recon[:, proj]
#         # print(np.shape(trng_dataseries_recon))

#         samples_uniform_recon = np.array([ecdf_func(xnew) for xnew in trng_dataseries_recon])

#         samples_uniform_target = samples_uniform[1:]
#         samples_uniform_recon = samples_uniform_recon[:-1]

#         # print(np.shape(samples_uniform_target))
#         # print(np.shape(samples_uniform_recon))

#         # percentage = np.sum(samples_uniform_target==samples_uniform_recon) / len(samples_uniform_target)
#         # print("percentage = {:}".format(percentage))

#         corr = np.corrcoef(samples_uniform_target,samples_uniform_recon)
#         print(corr)

#         data = {
#             "proj":proj,
#             "noise_level":noise_level,
#             "noise_var":noise_var,
#             "trng_data_noisy":trng_data_noisy,
#             "trng_data_recon":trng_data_recon,
#             "ecdf_func":ecdf_func,
#             "samples_uniform_target":samples_uniform_target,
#             "samples_uniform_recon":samples_uniform_recon,
#             "corr":corr,
#         }
#         data_path = saving_path_data + '/CrackingWithNoiseTest_proj{:}_lag{:}_noise_level{:}.pickle'.format(proj, trng_lag, noise_level)
#         utils.saveData(data, data_path, add_file_format=False)





# import matplotlib as mpl
# import matplotlib.pyplot as plt
# from matplotlib import rc
# rc('text', usetex=True)
# rc('font', family='serif')

# # for proj in range(D):
# for proj in [1]:
#     for noise_level in [0, 1, 5, 10]:
#         corr_values = []
#         for trng_lag in lagtimes:
#             data_path = saving_path_data + '/CrackingWithNoiseTest_proj{:}_lag{:}_noise_level{:}.pickle'.format(proj, trng_lag, noise_level)
#             data = utils.loadData(data_path, add_file_format=False)
#             corr = data["corr"][0,1]
#             corr_values.append(corr)

#         plt.plot(lagtimes, corr_values, label="Noise level {:}".format(noise_level))

# # plt.plot(lagtimes, Zcritical * np.ones_like(lagtimes), "r--")
# # plt.ylim([0, 2*Zcritical])

# figure_path = saving_path_figures + "/CrackingWithNoiseTest_proj{:}.pdf".format(proj)
# plt.legend(loc="upper left", bbox_to_anchor=(1.05, 1), borderaxespad=0.)
# plt.tight_layout()
# plt.savefig(figure_path)
# plt.close()










