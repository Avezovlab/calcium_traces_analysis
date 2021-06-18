#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 18:02:46 2021

@author: pierre

Plot individual traces contained in a given folder.
"""

from os.path import join, isfile, isdir
from os import listdir, makedirs
import pickle

from numpy import array, mean, max, std, histogram, correlate, argmax, arange, sqrt, cumsum
from numpy.polynomial.polynomial import polyfit

from matplotlib import pyplot as plt
from scipy.signal import find_peaks


from scipy.stats import pearsonr

from math import isnan

plt.ioff()


base_dir = "/mnt/data2/calcium_incucyte/Cas9_NMDA/traces"
out_dir = "/tmp/a"

if not isdir(out_dir):
    makedirs(out_dir)

DT = 0.33
pxsize = 0.2646 #mm

all_files =  [f for f in listdir(base_dir) if f.startswith("traces")]
for cpt,fname in enumerate(all_files):
    print("Processing[{}/{}]: {}".format(cpt, len(all_files), fname))
    traces = []
    with open(join(base_dir, fname), 'r') as f:
        for line in f:
            line = line.rstrip("\n").split(" ")
            traces.append([float(e) for e in line[1].split(",")])

    avg_sig = [mean([e[i] for e in traces]) for i in range(len(traces[0]))]

    plt.figure()
    for trace in traces:
        plt.plot(array(range(len(trace))) * DT, trace, linewidth=0.2)
    plt.plot(array(range(len(avg_sig))) * DT, avg_sig, 'k', linewidth=2)
    plt.ylim([0, 150])
    plt.ylabel("Intensity (AU)")
    plt.xlabel("Time (s)")
    plt.title("n = {} traces".format(len(traces)))
    plt.savefig(join(out_dir, fname[:-len("csv")] + ".png"), dpi=300)
    plt.close()
