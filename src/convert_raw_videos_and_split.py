from ij import IJ
import sys

import os
from os import path
sys.path.append("/mnt/data2/calcium_incucyte/PP_VD_Prop_191021/data")
from analysis_config import *


if not path.isdir(out_dir):
	os.makedirs(out_dir)

filenames = []
for root, dirs, files, in os.walk(in_dir):
	filenames.extend([f for f in files if f.endswith('.avi')])

for i, fname in enumerate(filenames):
	out_fname = fname.split("_")
	if out_fname[0] != "VID":
		out_fname[0] = "VID"
	out_fname = ("_".join(out_fname))[:-len(".avi")]

	print(path.join(root, fname))
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
