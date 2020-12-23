#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 16:15:58 2020

@author: pierre
"""

from skimage.measure import label, regionprops
from skimage.external.tifffile import imread, imsave
from numpy import zeros, array, mean, sum, sqrt, argmin
from os import listdir, makedirs
from os.path import join, isdir

pxsize = 0.2646 #mm
DT = 0.33 #s

dist_th = 50

cache_dir = "/mnt/data2/mosab_incucyte/processed/"


well_id = "B3_1"
timestamps = ["01d19h48m", "02d10h48m", "02d20h21m"]

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
    objs_cents[timestamp] = []
    with open(fname, 'r') as f:
        for line in f:
            line = line.strip("\n").split(" ")
            tmp = array([[float(e.split(",")[0]), float(e.split(",")[1])] for e in line[1:]])
            objs_cents[timestamp].append(mean(tmp, axis=0))

tab = []

seen = [[0] * len(traces[ts]) for ts in timestamps]
for k in range(len(timestamps) - 1):
    ts = timestamps[k]
    tss = timestamps[k+1]
    for l in range(len(objs_cents[ts])):
        if not seen[k][l]:
            seen[k][l] = 1
            tab.append([objs_cents[ts][l]])

        dists = [sqrt(sum((objs_cents[ts][l] - c)**2)) for c in objs_cents[tss]]
        min_d_idx = argmin(dists)

        if dists[min_d_idx] < dist_th:
            tab[-1].append(objs_cents[tss][min_d_idx])
            seen[k+1][min_d_idx] = 1