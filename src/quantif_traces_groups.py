#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  7 19:49:09 2021

@author: pierre
"""

import sys
from os.path import join, isdir
from os import listdir, makedirs

from numpy import mean, max, std, correlate, arange, argsort, isnan, array

from matplotlib import pyplot as plt
from scipy.signal import find_peaks
from scipy.interpolate import interp1d

from scipy.signal import savgol_filter
from scipy.stats import f_oneway

def pool_data_avg_per_well(d, group, timestamps, fill_missing):
    res = {}
    for k in d[group].keys():
        res[k] = [[], []]
        for l, ts in enumerate(timestamps):
            if ts in d[group][k] and len(d[group][k][ts]) > 0:
                res[k][0].append(l+1)
                res[k][1].append(d[group][k][ts])
            elif fill_missing:
                res[k][0].append(l+1)
                res[k][1].append(0)

    return res

def pool_data_timecat(d, timestamps, group, time_groups, fill_missing):
    ks = list(d[group].keys())
    res = {i:([0] * len(ks)) for i in set(time_groups)}
    cpts = {i:([0] * len(ks)) for i in set(time_groups)}
    for k in range(len(time_groups)):
        ts = timestamps[k]
        for j in range(len(ks)):
            vals = d[group][ks[j]]
            if ts in vals.keys() and len(vals[ts]) > 0:
                res[time_groups[k]][j] += mean(vals[ts])
                cpts[time_groups[k]][j] += 1
            elif fill_missing:
                res[time_groups[k]][j] += 0
                cpts[time_groups[k]][j] += 1

    for k in res.keys():
        for j in range(len(res[k])):
            res[k][j] /= cpts[k][j]

    return res

def pool_data(d, timestamps, groups, pool_time, fill_missing):
    xs = []
    dat = []
    if pool_times:
        for gk in groups.keys():
            xs.append(gk)
            dat.append([])
            for ts in timestamps:
                for vals in d[gk].values():
                    if ts in vals.keys() and len(vals[ts]) > 0:
                        dat[-1].extend(vals[ts])
                    elif fill_missing:
                        dat[-1].extend([0])
    else:
        for gk in groups.keys():
            dat.append([])
            for ts in timestamps:
                dat[-1].append([])
                if gk == "":
                    xs.append(ts)
                else:
                    xs.append(ts + "_" + gk)
                for vals in d[gk].values():
                    if ts in vals.keys() and len(vals[ts]) > 0:
                        dat[-1][-1].extend(vals[ts])
                    elif fill_missing:
                        dat[-1][-1].extend([0])
    return (xs, dat)

#plt.ioff()

pk_smooth_factor = 1.3

#base_dir = "/mnt/data2/calcium_incucyte/201228_Tasuku_MA24_16h cold shock/processed"
#out_dir = "/tmp/a/201228_Tasuku_MA24_16h cold shock"
#excluded_ts = []
#groups = {"WT": ["B2", "B3", "C2", "C3", "D2", "D3"],
#          "RTN3OE": ["B4", "C4", "D4"]}
#groups_col = {"WT": 'k', "RTN3OE": 'r'}
# time_groups = []
# time_groups_names = {}


# base_dir = "/mnt/data2/incucyte_calcium_cold_shock/MA9(24.8)_cold_shock2_19.9.20/processed"
# out_dir = "/tmp/a/MA9(24.8)_cold_shock2_19.9.20"
# excluded_ts = ["02d05h15m", "02d09h15m", "02d13h15m", "02d17h15m"]
# groups = {"": ["B2", "B3", "B4", "C2", "C3", "C4"]}
# time_groups = [1, 1, 1, 2, 2, 2, 2, 3, 3, 3]
# time_groups_names = {1:"before", 2:"schock", 3:"after"}
# groups_col = {"": 'k'}
# fill_missing = False
# plot_indiv_traces = False


# base_dir = "/mnt/data2/incucyte_calcium_cold_shock/MA12(18.9.20)_Cold shock3_12.10.20/processed"
# out_dir = "/tmp/a/MA12(18.9.20)_Cold shock3_12.10.20"
# excluded_ts = ["00d16h00m", "01d00h00m", "02d14h21m", "02d20h21m"]
# groups = {"": ["B3", "B4", "C2", "C3", "C4"]}
# groups_col = {"": 'k'}
# time_groups = [1, 1, 1, 2, 2 ,2, 2, 2, 3, 3, 3]
# time_groups_names = {1:"before", 2:"schock", 3:"after"}
# fill_missing = False
# plot_indiv_traces = True


# base_dir = "/mnt/data2/calcium_incucyte/MA11(31.8.20)_CPA_TG_WT _25.9.20/processed"
# excluded_ts = ["01d16h06m", "01d22h06m"]
# groups = {"CPA": ["B2"]}
# out_dir = "/tmp/a/MA11(31.8.20)_CPA_TG_WT_25.9.20/CPA"
# # groups = {"TG5µM": ["C3"]}
# # out_dir = "/tmp/a/MA11(31.8.20)_CPA_TG_WT_25.9.20/TG5"
# # groups = {"TG0.5µM": ["D3"]}
# # out_dir = "/tmp/a/MA11(31.8.20)_CPA_TG_WT_25.9.20/TG0p5"
# groups_col = {"CPA": 'k', "TG5µM": 'r', "TG0.5µM": "m"}
# time_groups = []
# time_groups_names = {}
# fill_missing = True
# plot_indiv_traces = True


# base_dir = "/mnt/data2/calcium_incucyte/Cas9_NMDA/processed"
# out_dir = "/tmp/a/Cas9_NMDA"
# excluded_ts = []
# groups = {"NMDA": ["B3"]}
# time_groups = []
# time_groups_names = {}
# fill_missing = True
# plot_indiv_traces = True


# base_dir = "/mnt/data2/calcium_incucyte/Cas9_AMPA/processed"
# out_dir = "/tmp/a/Cas9_AMPA"
# excluded_ts = []
# groups = {"AMPA": ["C3"]}
# time_groups = []
# time_groups_names = {}
# fill_missing = True
# plot_indiv_traces = True

# base_dir = "/mnt/nmve/vd_cyrcadyan/traces"
# out_dir = "/tmp/a/vd_cyrcadyan"
# excluded_ts = []
# groups = {"": ["C3", "C4", "C5"]}
# groups_col = {"": 'k'}
# time_groups = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
# time_groups_names = {1: "aa"}
# fill_missing = False
# plot_indiv_traces = False
# pk_smooth_factor=1.1

#these values must be present in the config file
proc_dir = None
pk_smooth_factor = None

sys.path.append("/mnt/data2/calcium_incucyte/ML/241121")
from analysis_config import *
excluded_ts = []
groups = {"chol" : ["C6", "C7", "C8"],
          "asyn" : ["D3", "D4", "D5"], "chol + asyn" : ["D6", "D7", "D8"], 
          "mbcd + asyn + chol" : ["E3", "E4", "E5"], "U188A + asyn + chol" : ["E6", "E7", "E8"],
          "U188A" : ["F3", "F4", "F5"], "cntrl" : ["F6", "F7", "F8"]}

#set1 matplotlib qualitative colormap
groups_col = {'chol': [228 / 255, 26 / 255, 28 / 255],
              'asyn': [55 / 255, 126 / 255, 184 / 255],
              'chol + asyn': [77 / 255, 175 / 255, 74 / 255],
              'mbcd + asyn + chol': [152 / 255, 78 / 255, 183 / 255],
              'U188A + asyn + chol': [255 / 255, 127 / 255, 0 / 255],
              'U188A': [166 / 255, 86 / 255, 40 / 255],
              'cntrl': [0, 0, 0]}
time_groups = [1] * 6

time_groups_names = {1: "aa"}
fill_missing = False
plot_indiv_traces = False

DT = 0.33
pxsize = 0.2646 #mm



gkeys = list(groups.keys())
pool_times = False
min_traces = 550 #50

timestamps = set()


base_dir = join(proc_dir, "traces")
out_dir = join(proc_dir, "res")

if not isdir(out_dir):
    makedirs(out_dir)

all_files =  [f for f in listdir(base_dir) if f.startswith("traces")]
all_traces_cnt = {k:{} for k in groups.keys()}
all_avg_pks_cnt = {k:{} for k in groups.keys()}
all_max_ccor = {k:{} for k in groups.keys()}
all_avg_pks_freqs = {k:{} for k in groups.keys()}
all_avg_pks_amps = {k:{} for k in groups.keys()}
all_avg_pks_fwhm = {k:{} for k in groups.keys()}
for cpt,fname in enumerate(all_files):
    print("Processing[{}/{}]: {}".format(cpt + 1, len(all_files), fname))
    ts = fname[:-len(".csv")].split("_")[-3]

    if ts in excluded_ts:
        continue

    well = "_".join(fname.split("_")[1:3])
    cat = None
    for k,v in groups.items():
        if any([e in fname for e in v]):
            cat = k
            break
    if cat is None:
        print("No category found for {}".format(fname))
        continue

    if well not in all_avg_pks_cnt[cat]:
        all_traces_cnt[cat][well] = {}
        all_max_ccor[cat][well] = {}
        all_avg_pks_cnt[cat][well] = {}
        all_avg_pks_freqs[cat][well] = {}
        all_avg_pks_amps[cat][well] = {}
        all_avg_pks_fwhm[cat][well] = {}

    traces = []
    with open(join(base_dir, fname), 'r') as f:
        for line in f:
            line = line.rstrip("\n").split(" ")
            traces.append([float(e) for e in line[1].split(",")])

    all_traces_cnt[cat][well][ts] = [len(traces)]
    if len(traces) < min_traces:
        all_avg_pks_cnt[cat][well][ts] = [0]
        all_max_ccor[cat][well][ts] = []
        all_avg_pks_freqs[cat][well][ts] = []
        all_avg_pks_amps[cat][well][ts] = []
        all_avg_pks_fwhm[cat][well][ts] = []
        continue

    timestamps.add(ts)

    avg_trace = [mean([e[i] for e in traces]) for i in range(len(traces[0]))]
    Nfilt = len(avg_trace)
    if Nfilt % 2 == 0:
        Nfilt -= 1
    smoothed_trace = savgol_filter(avg_trace, Nfilt, 3, mode="mirror")

    cur_pks = find_peaks(avg_trace)[0]
    Mtr = mean(avg_trace)
    SDtr = std(avg_trace)

    cur_pks = [pk for pk in cur_pks if ((smoothed_trace[pk] > 2 and avg_trace[pk] > pk_smooth_factor*smoothed_trace[pk]) or (smoothed_trace[pk] < 2 and avg_trace[pk] > 3*smoothed_trace[pk]))]
    all_avg_pks_cnt[cat][well][ts] = [len(cur_pks)]
    all_avg_pks_freqs[cat][well][ts] = [1 / ((cur_pks[i+1] - cur_pks[i]) * DT) for i in range(0, len(cur_pks) - 1) if (cur_pks[i+1] - cur_pks[i]) * DT < 50]
    all_avg_pks_amps[cat][well][ts] = [avg_trace[cur_pks[i]] for i in range(0, len(cur_pks))]

    all_max_ccor[cat][well][ts] = []
    for i in range(len(traces)):
        ccov = correlate(traces[i] - mean(traces[i]), avg_trace - mean(avg_trace), mode='full')
        ccor = ccov / (len(traces[i]) * std(traces[i]) * std(avg_trace))
        if not isnan(max(ccor)):
            all_max_ccor[cat][well][ts].append(max(ccor))

    I1s = []
    I2s = []
    yhs = []
    all_avg_pks_fwhm[cat][well][ts] = []
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
        all_avg_pks_fwhm[cat][well][ts].append((I2 - I1) * DT)
        I1s.append(I1)
        I2s.append(I2)
        yhs.append(y_h)

    if plot_indiv_traces:
        plt.figure(figsize=(6,6))
        for trace in traces:
            plt.plot(array(range(len(trace))) * DT, trace, linewidth=0.2)
        plt.plot(array(range(len(avg_trace))) * DT, avg_trace, 'k', linewidth=2)
        plt.plot([e * DT for e in cur_pks], [avg_trace[idx] for idx in cur_pks], '*r')
        plt.plot(array(range(len(smoothed_trace))) * DT, smoothed_trace, 'b', linewidth=2)
        for i in range(len(I1s)):
            plt.plot([I1s[i] * DT, I2s[i] * DT], [yhs[i], yhs[i]], 'g')
        plt.ylim([0, 175])
        plt.ylabel("Intensity (AU)")
        plt.xlabel("Time (s)")
        plt.title("n = {} traces".format(len(traces)))
        plt.savefig(join(out_dir, fname[:-len(".csv")] + ".png"), dpi=300)#, bbox_inches='tight',pad_inches=0)
        plt.close()


timestamps = sorted(timestamps,
                    key=lambda e: int(e.split("d")[0]) * 1440 + int(e.split("d")[1].split("h")[0]) * 60 + int(e.split("h")[1].split("m")[0]))
timestamp_disps = range(len(timestamps))
groups_keys = list(groups.keys())

xs, dat = pool_data(all_traces_cnt, timestamps, groups, pool_times, fill_missing)
plt.figure(figsize=(7,7))
for k in range(len(groups.keys())):
    plt.fill_between(range(1, len(timestamps)+1),
                      [mean(e) -std(e) for e in dat[k]],
                      [mean(e) + std(e) for e in dat[k]], alpha=0.3, color=groups_col[groups_keys[k]])
    plt.plot(range(1, len(timestamps)+1), [mean(e) for e in dat[k]], color=groups_col[groups_keys[k]])
plt.xticks(ticks=range(1, len(timestamps) + 1), labels=timestamp_disps, rotation=0)
plt.xlabel('Days after treatment')
plt.ylabel('Average number of traces')
plt.savefig(join(out_dir, "traces_count.png"), dpi=300, bbox_inches='tight',pad_inches=0)


xs, dat = pool_data(all_avg_pks_cnt, timestamps, groups, pool_times, fill_missing)
plt.figure(figsize=(7,7))
for k in range(len(groups.keys())):
    # plt.fill_between(range(1, len(timestamps)+1),
    #                  [mean(e) -std(e) for e in dat[k]],
    #                  [mean(e) + std(e) for e in dat[k]], alpha=0.3)
    plt.plot(range(1, len(timestamps)+1), [mean(e) for e in dat[k]], color=groups_col[groups_keys[k]])
plt.xticks(ticks=range(1, len(timestamps) + 1), labels=timestamp_disps, rotation=0)
plt.xlabel('Days after treatment')
plt.ylabel('Average number of peaks')
plt.savefig(join(out_dir, "peaks_count.png"), dpi=300, bbox_inches='tight',pad_inches=0)


xs, dat = pool_data(all_avg_pks_freqs, timestamps, groups, pool_times, fill_missing)
plt.figure(figsize=(7,7))
for k in range(len(groups.keys())):
    # plt.fill_between(range(1, len(timestamps)+1),
    #                   [mean(e) -std(e) for e in dat[k]],
    #                   [mean(e) + std(e) for e in dat[k]], alpha=0.3)
    plt.plot(range(1, len(timestamps)+1), [mean(e) for e in dat[k]], color=groups_col[groups_keys[k]])
plt.xticks(ticks=range(1, len(timestamps) + 1), labels=timestamp_disps, rotation=0)
plt.xlabel('Days after treatment')
plt.ylabel('Peak Frequency (Hz)')
plt.savefig(join(out_dir, "peaks_freq.png"), dpi=300, bbox_inches='tight',pad_inches=0)

xs, dat = pool_data(all_avg_pks_amps, timestamps, groups, pool_times, fill_missing)
plt.figure(figsize=(7,7))
for k in range(len(groups.keys())):
    # plt.fill_between(range(1, len(timestamps)+1),
    #                   [mean(e) -std(e) for e in dat[k]],
    #                   [mean(e) + std(e) for e in dat[k]], alpha=0.3)
    plt.plot(range(1, len(timestamps)+1), [mean(e) for e in dat[k]], color=groups_col[groups_keys[k]])
plt.xticks(ticks=range(1, len(timestamps) + 1), labels=timestamp_disps, rotation=0)
plt.xlabel('Days after treatment')
plt.ylabel('Peak amplitude (AU)')
plt.savefig(join(out_dir, "peaks_amplitude.png"), dpi=300, bbox_inches='tight',pad_inches=0)


xs, dat = pool_data(all_avg_pks_fwhm, timestamps, groups, pool_times, fill_missing)
plt.figure(figsize=(7,7))
for k in range(len(groups.keys())):
    plt.fill_between(range(1, len(timestamps)+1),
                      [mean(e) -std(e) for e in dat[k]],
                      [mean(e) + std(e) for e in dat[k]], alpha=0.3, color=groups_col[groups_keys[k]])
    plt.plot(range(1, len(timestamps)+1), [mean(e) for e in dat[k]], color=groups_col[groups_keys[k]])
plt.xticks(ticks=range(1, len(timestamps) + 1), labels=timestamp_disps, rotation=0)
plt.xlabel('Days after treatment')
plt.ylabel('Peak width (s)')
plt.savefig(join(out_dir, "peaks_width.png"), dpi=300, bbox_inches='tight',pad_inches=0)

plt.figure(figsize=(7,7))
for i, k in enumerate(groups_keys):
    vs = []
    for w,v in all_avg_pks_fwhm[k].items():
        vs.append([mean(v[ts]) for ts in timestamps])
    plt.plot(range(1, len(timestamps) + 1), [mean([e[i] for e in vs]) for i in range(6)], color=groups_col[groups_keys[i]])
    plt.fill_between(range(1, len(timestamps) + 1),
                     [mean([e[i] for e in vs]) - std([e[i] for e in vs]) for i in range(6)],
                     [mean([e[i] for e in vs]) + std([e[i] for e in vs]) for i in range(6)], alpha=0.3, color=groups_col[groups_keys[i]])
plt.xticks(ticks=range(1, len(timestamps) + 1), labels=timestamp_disps, rotation=0)
plt.xlabel('Days after treatment')
plt.ylabel('Average peak width (s)')
plt.savefig(join(out_dir, "peaks_width_notpooled.png"), dpi=300, bbox_inches='tight',pad_inches=0)



xs, dat = pool_data(all_max_ccor, timestamps, groups, pool_times, fill_missing)
plt.figure(figsize=(7,7))
for k in range(len(groups.keys())):
    # plt.fill_between(range(1, len(timestamps)+1),
    #                   [mean(e) -std(e) for e in dat[k]],
    #                   [mean(e) + std(e) for e in dat[k]], alpha=0.3)
    plt.plot(range(1, len(timestamps)+1), [mean(e) for e in dat[k]], color=groups_col[groups_keys[k]])
plt.xticks(ticks=range(1, len(timestamps) + 1), labels=timestamp_disps, rotation=0)
plt.xlabel('Days after treatment')
plt.ylabel('AVG max correlation')
plt.savefig(join(out_dir, "corr.png"), dpi=300, bbox_inches='tight',pad_inches=0)

plt.figure(figsize=(7,7))
for i, k in enumerate(groups_keys):
    vs = []
    for w,v in all_max_ccor[k].items():
        vs.append([mean(v[ts]) for ts in timestamps])
    plt.plot(range(1, len(timestamps)+1), [mean([e[i] for e in vs]) for i in range(6)], color=groups_col[groups_keys[i]])
    plt.fill_between(range(1, len(timestamps)+1),
                     [mean([e[i] for e in vs]) - std([e[i] for e in vs]) for i in range(6)],
                     [mean([e[i] for e in vs]) + std([e[i] for e in vs]) for i in range(6)], alpha=0.3, color=groups_col[groups_keys[i]])
plt.xticks(ticks=range(1, len(timestamps) + 1), labels=timestamp_disps, rotation=0)
plt.xlabel('Days after treatment')
plt.ylabel('AVG delta max correlation')
plt.savefig(join(out_dir, "corr_notpooled.png"), dpi=300, bbox_inches='tight',pad_inches=0)

plt.figure(figsize=(7,7))
for i, k in enumerate(groups_keys):
    vs = []
    for w,v in all_max_ccor[k].items():
        n = mean(v[timestamps[0]])
        vs.append([mean(v[ts] - n) for ts in timestamps])
    plt.plot(range(1, len(timestamps)+1), [mean([e[i] for e in vs]) for i in range(6)], color=groups_col[groups_keys[i]])
    plt.fill_between(range(1, len(timestamps)+1),
                     [mean([e[i] for e in vs]) - std([e[i] for e in vs]) for i in range(6)],
                     [mean([e[i] for e in vs]) + std([e[i] for e in vs]) for i in range(6)], alpha=0.3, color=groups_col[groups_keys[i]])
plt.xticks(ticks=range(1, len(timestamps) + 1), labels=timestamp_disps, rotation=0)
plt.xlabel('Days after treatment')
plt.ylabel('AVG delta max correlation')
plt.savefig(join(out_dir, "corr_delta_notpooled.png"), dpi=300, bbox_inches='tight',pad_inches=0)



