#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 16:46:04 2021

@author: pierre
"""

from os.path import join, isdir
from os import listdir, makedirs

from numpy import array, mean, max, std, arange, argsort
from numpy.random import rand

from matplotlib import pyplot as plt
from scipy.signal import find_peaks
from scipy.interpolate import interp1d

from scipy.signal import savgol_filter

#plt.ioff()


#base_dir = "/mnt/data2/incucyte_calcium_cold_shock/MA12(18.9.20)_Cold shock3_12.10.20/processed"
#out_dir = "/tmp/a/MA12(18.9.20)_Cold shock3_12.10.20"

#base_dir = "/mnt/data2/incucyte_calcium_cold_shock/MA9(24.8)_cold_shock2_19.9.20/processed"
#out_dir = "/tmp/a/MA9(24.8)_cold_shock2_19.9.20"

import sys
sys.path.append("/mnt/data2/calcium_incucyte/PP_VD_Prop_191021/data")
from analysis_config import *

DT = 0.33
pxsize = 0.2646 #mm

plot_indiv_traces = True

base_dir = join(out_dir, "traces")
out_dir = join(out_dir, "imgs")
if not isdir(out_dir):
    makedirs(out_dir)


all_files =  [f for f in listdir(base_dir) if f.startswith("traces")]
timestamps = set()
pks_cnt = {}
all_avg_pks_freqs = {}
all_avg_pks_amps = {}
all_avg_pks_fwhm = {}
for cpt,fname in enumerate(all_files):
    print("Processing[{}/{}]: {}".format(cpt + 1, len(all_files), fname))
    ts = fname[:-len(".csv")].split("_")[3]
    timestamps.add(ts)

    well = "_".join(fname.split("_")[1:3])
    if well not in pks_cnt:
        pks_cnt[well] = {}
        all_avg_pks_freqs[well] = {}
        all_avg_pks_amps[well] = {}
        all_avg_pks_fwhm[well] = {}

    if ts not in pks_cnt[well]:
        pks_cnt[well][ts] = []
        all_avg_pks_freqs[well][ts] = []
        all_avg_pks_amps[well][ts] = []
        all_avg_pks_fwhm[well][ts] = []

    traces = []
    with open(join(base_dir, fname), 'r') as f:
        for line in f:
            line = line.rstrip("\n").split(" ")#?
            traces.append([float(e) for e in line[1].split(",")])
    avg_trace = [mean([e[i] for e in traces]) for i in range(len(traces[0]))]
    smoothed_trace = savgol_filter(avg_trace, 179, 3, mode='mirror')

    cur_pks = find_peaks(avg_trace)[0]
    Mtr = mean(avg_trace)
    SDtr = std(avg_trace)

    cur_pks = [pk for pk in cur_pks if avg_trace[pk] > 1.3*smoothed_trace[pk]]
    
    pks_cnt[well][ts] = len(cur_pks)
    all_avg_pks_freqs[well][ts] = [1 / ((cur_pks[i+1] - cur_pks[i]) * DT) for i in range(0, len(cur_pks) - 1) if (cur_pks[i+1] - cur_pks[i]) * DT < 50]
    all_avg_pks_amps[well][ts] = [avg_trace[cur_pks[i]] for i in range(0, len(cur_pks))]

    I1s = []
    I2s = []
    yhs = []
    all_avg_pks_fwhm[ts] = []
    prev_pk = 0
    next_pk = 1
    for cpt_pk, pk in enumerate(cur_pks):
        if cpt_pk > 0:
            prev_pk = cur_pks[cpt_pk - 1]
        if cpt_pk < len(cur_pks) - 1:
            next_pk = cur_pks[cpt_pk+1]
        else:
            next_pk = max(len(avg_trace))

        prev_min = argsort(avg_trace[prev_pk:pk])[0] + prev_pk
        next_min = argsort(avg_trace[pk:next_pk])[0] + pk

        y = avg_trace[pk]
        y_h = y - (y - max([avg_trace[prev_min], avg_trace[next_min]])) / 2
        f = interp1d(range(len(avg_trace)), avg_trace)
        xs = arange(max([pk - 10, 0]), min([pk + 10, len(avg_trace) - 1]), 0.1)
        I = argsort(abs(f(xs) - y_h))
        I1 = [xs[e] for e in I if xs[e] < pk and xs[e] > prev_min][0]
        I2 = [xs[e] for e in I if xs[e] > pk and xs[e] < next_min][0]
        all_avg_pks_fwhm[well][ts].append((I2 - I1) * DT)
        I1s.append(I1)
        I2s.append(I2)
        yhs.append(y_h)

    if plot_indiv_traces:
        plt.figure(figsize=(6,2))
        for trace in traces:
            plt.plot(array(range(len(trace))) * DT, trace, linewidth=0.2)
        plt.plot(array(range(len(avg_trace))) * DT, avg_trace, 'k', linewidth=2)
        plt.plot([e * DT for e in cur_pks], [avg_trace[idx] for idx in cur_pks], '*r')
        #plt.plot(array(range(len(smoothed_trace))) * DT, smoothed_trace, 'b', linewidth=2)
        for i in range(len(I1s)):
            plt.plot([I1s[i] * DT, I2s[i] * DT], [yhs[i], yhs[i]], 'g')
        plt.ylim([0, 175])
        plt.ylabel("Intensity (AU)")
        plt.xlabel("Time (s)")
        plt.title("n = {} traces".format(len(traces)))
        #plt.savefig(join(out_dir, fname[:-len(".csv")] + ".png"), dpi=300)#, bbox_inches='tight',pad_inches=0)
        plt.savefig(join(out_dir, fname[:-len(".csv")] + ".svg"))
        plt.close()
        #assert(False)


timestamps = sorted(timestamps,
                    key=lambda e: int(e.split("d")[0]) * 1440 + int(e.split("d")[1].split("h")[0]) * 60 + int(e.split("h")[1].split("m")[0]))


min_v = min([min(e.values()) for e in pks_cnt.values()])
max_v = max([max(list(e.values())) for e in pks_cnt.values()])
for exp in pks_cnt.keys():
    plt.figure(figsize=(6,6))
    plt.plot(range(1, len(timestamps) + 1), [pks_cnt[exp][ts] for ts in timestamps])
    plt.xticks(ticks=range(1, len(timestamps) + 1), labels=timestamps, rotation=45)
    plt.ylim([min_v, max_v])
    plt.ylabel('Number of peaks')
    plt.savefig(join(out_dir, exp + "_peaks_count.png"), dpi=300)#, bbox_inches='tight',pad_inches=0)


plt.figure(figsize=(6,6))
for exp in pks_cnt.keys():
    plt.plot(range(1, len(timestamps) + 1), [pks_cnt[exp][ts] for ts in timestamps])
plt.xticks(ticks=range(1, len(timestamps) + 1), labels=timestamps, rotation=45)
plt.ylabel('Number of peaks')
plt.legend(pks_cnt.keys())
#plt.savefig(join(out_dir, "peaks_count.png"), dpi=300)#, bbox_inches='tight',pad_inches=0)
plt.savefig(join(out_dir, "peaks_count.pdf"))#, bbox_inches='tight',pad_inches=0)

plt.figure(figsize=(6,6))
p = plt.boxplot([all_avg_pks_freqs[ts] for ts in timestamps])
for i,ts in enumerate(timestamps):
    plt.plot(i + 1 + (rand(len(all_avg_pks_freqs[ts])) - 0.5) * 0.2, all_avg_pks_freqs[ts], 'xk')
plt.xticks(ticks=range(1, len(timestamps) + 1), labels=timestamps, rotation=45)
plt.ylabel('Peak Frequency (Hz)')
plt.savefig(join(out_dir, "peaks_freq.png"), dpi=300)#, bbox_inches='tight',pad_inches=0)

plt.figure(figsize=(6,6))
p = plt.boxplot([all_avg_pks_amps[ts] for ts in timestamps])
for i,ts in enumerate(timestamps):
    plt.plot(i + 1 + (rand(len(all_avg_pks_amps[ts])) - 0.5) * 0.2, all_avg_pks_amps[ts], 'xk')
plt.xticks(ticks=range(1, len(timestamps) + 1), labels=timestamps, rotation=45)
plt.ylabel('Peak amplitude (AU)')
plt.savefig(join(out_dir, "peaks_amp.png"), dpi=300)#, bbox_inches='tight',pad_inches=0)

plt.figure(figsize=(6,6))
p = plt.boxplot([all_avg_pks_fwhm[ts] for ts in timestamps])
for i,ts in enumerate(timestamps):
    plt.plot(i + 1 + (rand(len(all_avg_pks_fwhm[ts])) - 0.5) * 0.2, all_avg_pks_fwhm[ts], 'xk')
plt.xticks(ticks=range(1, len(timestamps) + 1), labels=timestamps, rotation=45)
plt.ylabel('Full width at half maximum (s)')
plt.savefig(join(out_dir, "peaks_width.png"), dpi=300)#, bbox_inches='tight',pad_inches=0)