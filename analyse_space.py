#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 15:39:44 2020

@author: pierre
"""

from os.path import join, isfile
import pickle

from numpy import array, loadtxt, mean, max, std, histogram, correlate, argmax, arange, sqrt, cumsum, histogram, zeros

from matplotlib import pyplot as plt
from scipy.signal import find_peaks

from skimage.external.tifffile import imsave


from scipy.stats import pearsonr

from math import isnan

cache_dir = "/mnt/data2/mosab_incucyte/processed/"

DT = 0.33
pxsize = 0.2646 #mm

well_id = "B3_1"
timestamps = ["01d19h48m", "02d10h48m", "02d20h21m"]

colors = {timestamps[0]:"#00FF49FF", timestamps[1]: "#FF4900FF", timestamps[2]: "#4900FFFF"}

traces = {}
objs_cents = {}
for timestamp in timestamps:
    fname = "/mnt/data2/mosab_incucyte/processed/traces_" + well_id + "_" + timestamp + ".csv"
    traces[timestamp] = []
    with open(fname, 'r') as f:
        for line in f:
            line = line.rstrip("\n").split(" ")
            traces[timestamp].append([float(e) for e in line[1].split(",")])

    fname = "/mnt/data2/mosab_incucyte/processed/objects_" + well_id + "_" + timestamp + ".csv"

    obj_img = zeros((1152, 1536))
    objs_cents[timestamp] = []
    with open(fname, 'r') as f:
        for line in f:
            line = line.strip("\n").split(" ")
            tmp = array([[float(e.split(",")[0]), float(e.split(",")[1])] for e in line[1:]])
            obj_img[tmp[:,0].astype("int"), tmp[:,1].astype("int")] = int(line[0])
            objs_cents[timestamp].append(mean(tmp, axis=0))
    imsave("/tmp/objects_labeled_" + well_id + "_" + timestamp + ".tif", obj_img)
    
