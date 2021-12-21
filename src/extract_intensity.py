#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 31 01:27:10 2020

@author: pierre

Excract individual cell's calcium traces from videos.

Inputs:
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

import sys
sys.path.append("/mnt/data2/calcium_incucyte/PP_VD_Prop_191021/data")
from analysis_config import *



vid_fmt = "VID_{}_{}_{}.tif"
mask_fmt = "MASK_{}_{}.tif"
out_obj_fmt = "objects_{}_{}.csv"
out_fmt = "{}_{}_{}_{}.csv"

files = listdir(join(in_dir, "tiffs"))
for cpt,fname in enumerate(files):
    timestamp = fname.split("_")[3]
    well_id = "_".join(fname.split("_")[1:3])
    substack="_".join(fname.split("_")[4:])[:-len(".tif")]

    vid_fname = join(in_dir, "tiffs", vid_fmt.format(well_id, timestamp, substack))
    mask_fname = join(in_dir, "masks", mask_fmt.format(well_id, timestamp))
    fr_mask_fname = join(in_dir, "fr_masks", mask_fmt.format(well_id, timestamp))


    if not isfile(vid_fname) or not isfile(mask_fname):
        print("Missing file [{}/{}]: {}".format(cpt+1, len(files), fname))
        continue


    out_obj_fname = join(in_dir, "traces", out_obj_fmt.format(well_id, timestamp))
    out_trace_fname = join(in_dir, "traces", out_fmt.format("traces", well_id, timestamp, substack))

    if isfile(out_trace_fname):
        print("Skipped[{}/{}]: {}".format(cpt+1, len(files), fname))
        continue
    print("Processing[{}/{}]: {}".format(cpt+1, len(files), fname))

    stck = mean(imread(vid_fname), axis=3)#mean color, to not differentiate among RGB
    mask = imread(mask_fname)

    if mask.ndim == 3:
        mask = mask[:,:,0] > 0#?

    labs = label(mask)
    labs_props = regionprops(labs)

    cell_Is = []
    for lab_prop in labs_props:
        cell_Is.append([])#list of lists
        for k in range(stck.shape[0]):#add to the list just created
            cell_Is[-1].append(mean(stck[k, lab_prop.coords[:,0], lab_prop.coords[:,1]]))#mean value of all pixel of an object

    with open(out_obj_fname, 'w') as f:
        for i,lab_prop in enumerate(labs_props):
            f.write(str(i) + " " + " ".join([str(lab_prop.coords[i,0]) + "," + str(lab_prop.coords[i,1]) for i in range(lab_prop.coords.shape[0])]) + "\n")#should obj and traces be the same for the two halves?? obj no, trace yes

    with open(out_trace_fname, 'w') as f:
        for i,trace in enumerate(cell_Is):
            f.write(str(i) + " " + ",".join(["{:.3f}".format(e) for e in trace]) + "\n")#f=Floating point decimal format. 3= round to 3 places after the decimal point

    if isfile(fr_mask_fname):
        print("Processing FR mask [{}/{}]: {}".format(cpt+1, len(files), fname))
        out_fr_obj_fname = join(in_dir, "traces", out_fmt.format("fr_objects", well_id, timestamp, substack))
        out_fr_trace_fname = join(in_dir, "traces", out_fmt.format("fr_traces", well_id, timestamp, substack))

        fr_mask = imread(mask_fname)
        if fr_mask.ndim == 3:
            fr_mask = fr_mask[:,:,0] > 0
        fr_mask = fr_mask and mask#where they're both true

        fr_labs = label(fr_mask)
        fr_labs_props = regionprops(fr_labs)

        fr_cell_Is = []
        for fr_lab_prop in fr_labs_props:
            fr_cell_Is.append([])
            for k in range(stck.shape[0]):
                fr_cell_Is[-1].append(mean(stck[k, fr_lab_prop.coords[:,0], fr_lab_prop.coords[:,1]]))

        with open(out_fr_obj_fname, 'w') as f:
            for i, fr_lab_prop in enumerate(labs_props):
                f.write(str(i) + " " + " ".join([str(fr_lab_prop.coords[i,0]) + "," + str(fr_lab_prop.coords[i,1]) for i in range(fr_lab_prop.coords.shape[0])]) + "\n")

        with open(out_fr_trace_fname, 'w') as f:
            for i, fr_trace in enumerate(fr_cell_Is):
                f.write(str(i) + " " + ",".join(["{:.3f}".format(e) for e in fr_trace]) + "\n")


