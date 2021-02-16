#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 31 01:27:10 2020

@author: pierre

Excract individual cell's calcium traces from videos.

Inputs:
  base dir must conain two sub-direcies: videos and masks
  video files must be tif and start with "VID_"
  mask files must be black and white tif images with names starting with "MASK_"

Outputs: the script generates two outputs per video as CSV files:
  * objects_XX.csv that contains for each line the id and pixels corresponding to an "active object"
  * traces_XX.csv that contains on each line the ide and mean intensity values over each pixels (as defined in objects_XX.csv) for each frame of the video
"""

from skimage.measure import label, regionprops
from skimage.io import imread
from numpy import mean
from os.path import join, isfile
from os import listdir


base_dir = "XX"
out_dir = "XX"



files = listdir(join(base_dir, "videos"))
for cpt,fname in enumerate(files):
    timestamp = fname.split("_")[3][:-len(".tif")]
    well_id = "_".join(fname.split("_")[1:3])

    vid_fname = join(base_dir, "videos", "VID_" + well_id + "_" + timestamp + ".tif")
    mask_fname = join(base_dir, "masks", "MASK_" + well_id + "_" + timestamp + ".tif")

    out_obj_fname = join(out_dir, "objects_" + well_id + "_" + timestamp + ".csv")
    out_trace_fname = join(out_dir, "traces_" + well_id + "_" + timestamp + ".csv")

    if isfile(out_obj_fname) or isfile(out_trace_fname):
        print("Skipped[{}/{}]: {}".format(cpt, len(files), fname))
        continue
    print("Processing[{}/{}]: {}".format(cpt, len(files), fname))

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

    with open(out_obj_fname, 'w') as f:
        for i,lab_prop in enumerate(labs_props):
            f.write(str(i) + " " + " ".join([str(lab_prop.coords[i,0]) + "," + str(lab_prop.coords[i,1]) for i in range(lab_prop.coords.shape[0])]) + "\n")

    with open(out_trace_fname, 'w') as f:
        for i,trace in enumerate(cell_Is):
            f.write(str(i) + " " + ",".join(["{:.3f}".format(e) for e in trace]) + "\n")


