#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 12:21:15 2021

@author: pierre
"""

from skimage.io import imread
from skimage.measure import label
from skimage.morphology import binary_erosion

from skimage.color import label2rgb
import matplotlib.pyplot as plt

from random import randint

im = imread("/mnt/data2/calcium_incucyte/VD_old_neurons_111221/MASK_VID958_D6_1_09d00h00m_Simple Segmentation.tif")

im[im == 1] = 0
im[im == 2] = 1

im = binary_erosion(im)

labs = label(im)


plt.figure(figsize=(20,20))
plt.imshow(label2rgb(labs, bg_label=0))
plt.savefig("/tmp/labs.png")N