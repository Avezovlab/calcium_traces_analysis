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


in_dir = "E:/vd_cyrcadyan"
out_dir = "E:/out/videos"

Nsplit = 2

filenames = []
for root, dirs, files, in os.walk(in_dir):
	filenames.extend([f for f in files if f.endswith('.avi')])

for i, fname in enumerate(filenames):
	out_fname = fname.split("_")
	if out_fname[0] != "VID":
		out_fname[0] = "VID"
	out_fname = ("_".join(out_fname))[:-len(".avi")]
	
	IJ.open(path.join(root, fname))
	imp = IJ.getImage()

	Nslices = imp.getNSlices()

	out_fmt = "{}_{}_{}.tif"
	if Nsplit == 1:
		out_path = path.join(out_dir, out_fmt.format(out_fname, 1, Nslices))
		IJ.save(out_path)
		imp.close()
	else:
		delta = Nslices / Nsplit
		for n in range(Nsplit):
			start = n * delta + 1
			end = (n + 1) * delta
			IJ.run("Make Substack...", "slices=" + str(start) + "-" + str(end))
			out_path = path.join(out_dir, out_fmt.format(out_fname, start, end))
			if not path.isfile(out_path):
				print("Processing[{}/{}:{}]: {}".format(i+1, len(filenames), n, fname))
				IJ.save(out_path)
			else:
				print("Skipped[{}/{}:{}]: {}".format(i+1, len(filenames), n, fname))
			imp = IJ.getImage()
			imp.close()
	imp = IJ.getImage()
	imp.close()
print("done")