#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 31 01:27:10 2020

@author: pierre
"""

from skimage.measure import label, regionprops
from skimage.external.tifffile import imread, imsave
from numpy import zeros, array, mean
from os import listdir, makedirs
from os.path import join, isdir

DT = 0.33

base_dir = "/mnt/data2/mosab_incucyte/data/"
out_dir = "/mnt/data2/mosab_incucyte/processed/"

well_id = "B3_1"
timestamps = ["01d19h48m", "02d10h48m", "02d20h21m"]

for timestamp in timestamps:
    vid_fname = join(base_dir, "videos", "VID_" + well_id + "_" + timestamp + ".tif")
    mask_fname = join(base_dir, "masks", "MASK_" + well_id + "_" + timestamp + ".tif")

    stck = mean(imread(vid_fname), axis=3)
    mask = imread(mask_fname)

    mask = mask[:,:,0] > 0

    labs = label(mask)
    labs_props = regionprops(labs)

    cell_Is = []
    for lab_prop in labs_props:
        cell_Is.append([])
        for k in range(stck.shape[0]):
            cell_Is[-1].append(mean(stck[k, lab_prop.coords[:,0], lab_prop.coords[:,1]]))

    with open(join(out_dir, "objects_" + well_id + "_" + timestamp + ".csv"), 'w') as f:
        for i,lab_prop in enumerate(labs_props):
            f.write(str(i) + " " + " ".join([str(lab_prop.coords[i,0]) + "," + str(lab_prop.coords[i,1]) for i in range(lab_prop.coords.shape[0])]) + "\n")

    with open(join(out_dir, "traces_" + well_id + "_" + timestamp + ".csv"), 'w') as f:
        for i,trace in enumerate(cell_Is):
            f.write(str(i) + " " + ",".join(["{:.3f}".format(e) for e in trace]) + "\n")


