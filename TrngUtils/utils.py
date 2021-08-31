#!/usr/bin/env python
# # -*- coding: utf-8 -*-

"""Created by: Vlachas Pantelis, CSE-lab, ETH Zurich
"""
#!/usr/bin/env python

import pickle
import numpy as np

def saveData(data, data_path, add_file_format=False):
    if add_file_format: data_path += ".pickle"
    with open(data_path, "wb") as file:
        # Pickle the "data" dictionary using the highest protocol available.
        pickle.dump(data, file, pickle.HIGHEST_PROTOCOL)
        del data
    return 0

def loadData(data_path, add_file_format=False):
    if add_file_format: data_path += ".pickle"
    try:
        with open(data_path, "rb") as file:
            data = pickle.load(file)
    except Exception as inst:
        print("[utils] Datafile\n {:s}\nNOT FOUND.".format(data_path))
        raise ValueError(inst)
    return data

def getNumberOfBins(N, rule="rice"):
    ######################
    ## Sturges rule
    ######################
    if rule == "sturges":
        nbins = int(1 + np.log2(N))
    ######################
    ## Rice rule
    ######################
    if rule == "rice":
        nbins = 2.0 * np.power(N, 1.0 / 3.0)
        nbins = int(np.ceil(nbins))
        nbins = int(np.max([2, nbins]))
    return nbins

def writeDataToByteArray(data, filename):
    """ Assert that the array is binary """
    assert np.all(np.logical_or(data==0, data==1))

    data = [int(d) for d in data]
    data = np.array(data)

    print("# ZEROS = {:}".format(np.sum(data==0)))
    print("# ONES = {:}".format(np.sum(data==1)))

    s = "".join(["{:}".format(d) for d in data])

    i = 0
    buffer = bytearray()
    while i < len(s):
        buffer.append( int(s[i:i+8], 2) )
        i += 8

    # now write your buffer to a file
    with open(filename, 'bw') as f:
        f.write(buffer)

    return 0



