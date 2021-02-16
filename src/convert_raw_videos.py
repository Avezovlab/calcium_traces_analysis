from fiji.plugin.trackmate import Model
from fiji.plugin.trackmate import Settings
from fiji.plugin.trackmate import TrackMate
from fiji.plugin.trackmate import SelectionModel
from fiji.plugin.trackmate import Logger
from fiji.plugin.trackmate.features.spot import SpotContrastAndSNRAnalyzerFactory, SpotIntensityAnalyzerFactory
from fiji.plugin.trackmate.detection import LogDetectorFactory
from fiji.plugin.trackmate.tracking.sparselap import SimpleSparseLAPTrackerFactory
from fiji.plugin.trackmate.tracking import LAPUtils
from ij import IJ, WindowManager
from fiji.plugin.trackmate import Spot
import fiji.plugin.trackmate.visualization.hyperstack.HyperStackDisplayer as HyperStackDisplayer
import fiji.plugin.trackmate.features.FeatureFilter as FeatureFilter
import sys
import fiji.plugin.trackmate.features.track.TrackDurationAnalyzer as TrackDurationAnalyzer
from ij.plugin import Commands
from ij.measure import Calibration
from ij.gui import Roi

import os
from os import path


in_dir = "/media/pierre/pierreDisk/MA12(18.9.20)_Cold shock3_12.10.20"
out_dir = "/mnt/data2/incucyte_calcium_cold_shock/MA12(18.9.20)_Cold shock3_12.10.20/data/videos"

force = False


filenames = []
for root, dirs, files, in os.walk(in_dir):
	filenames.extend([f for f in files if f.endswith('.avi')])

for i, fname in enumerate(filenames):
	out_fname = fname.split("_")
	if out_fname[0] != "VID":
		out_fname[0] = "VID"
	out_fname = "_".join(out_fname)
	
	out_path = path.join(out_dir, out_fname[:-len(".avi")] + ".tif")

	print(out_path)
	if not path.isfile(out_path):
		print("Processing[{}/{}]: {}".format(i, len(filenames), fname))
		IJ.open(path.join(root, fname))
		IJ.save(out_path)
		imp = IJ.getImage()
		imp.close()
	else:
		print("Skipped[{}/{}]: {}".format(i, len(filenames), fname))
print("done")