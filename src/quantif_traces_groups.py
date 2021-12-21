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
        for ts in timestamps:
            dat.append([])
            for gk in groups.keys():
                if gk == "":
                    xs.append(ts)
                else:
                    xs.append(ts + "_" + gk)
                for vals in d[gk].values():
                    if ts in vals.keys() and len(vals[ts]) > 0:
                        dat[-1].extend(vals[ts])
                    elif fill_missing:
                        dat[-1].extend([0])
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


sys.path.append("/mnt/data2/calcium_incucyte/PP_VD_Prop_191021/data")
from analysis_config import *

excluded_ts = []
groups = {"1": ["D6"], "2": ["D7"]}#, "CTRL2": ["E6"], "0.5_1": ["D7"], "0.5_2": ["D8"], "1.5_1": ["E7"], "1.5_2": ["E8"]}
groups_col = {"": 'k'}
time_groups = [1] * 15
time_groups_names = {1: "aa"}
fill_missing = False
plot_indiv_traces = False
pk_smooth_factor=1.1

if not isdir(out_dir):
    makedirs(out_dir)

DT = 0.33
pxsize = 0.2646 #mm



gkeys = list(groups.keys())
pool_times = False
min_traces = 50

timestamps = set()


base_dir = join(out_dir, "traces")
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
        continue

    if well not in all_avg_pks_cnt[cat]:
        all_traces_cnt[cat][well] = {}
        all_max_ccor[cat][well] = {}
        all_avg_pks_cnt[cat][well] = {}
        all_avg_pks_freqs[cat][well] ={}
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


# xs, dat = pool_data(all_traces_cnt, timestamps, groups, pool_times)
# v_x = []
# v_y = []
# for i in range(len(dat)):
#     v_x.extend([i] * len(dat[i]))
#     v_y.extend(dat[i])
# plt.figure(figsize=(7,7))
# p = sns.violinplot(x=v_x, y=v_y, cut=0)
# plt.ylabel('Number of traces')
# p.set_xticklabels(timestamps, rotation=90)
# plt.savefig(join(out_dir, "traces_count.png"), dpi=300, bbox_inches='tight',pad_inches=0)

#plt.figure(figsize=(7,7))
#p = plt.boxplot(dat)
#for i in range(len(dat)):
#    plt.plot(i + 1 + (rand(len(dat[i])) - 0.5) * 0.2, dat[i], 'x' + groups_col[gkeys[i % len(groups.keys())]])
#    plt.xticks(ticks=range(1, len(dat) + 1), labels=xs, rotation=90)
#plt.ylabel('Number of traces')
#plt.savefig(join(out_dir, "traces_count.png"), dpi=300, bbox_inches='tight',pad_inches=0)

# xs, dat = pool_data(all_avg_pks_cnt, timestamps, groups, pool_times)
# plt.figure(figsize=(7,7))
# p = plt.boxplot(dat)
# for i in range(len(dat)):
#     plt.plot(i + 1 + (rand(len(dat[i])) - 0.5) * 0.2, dat[i], 'x' + groups_col[gkeys[i % len(groups.keys())]])
# plt.xticks(ticks=range(1, len(dat) + 1), labels=xs, rotation=90)
# plt.ylabel('Number of peaks')
# plt.savefig(join(out_dir, "peaks_count.png"), dpi=300, bbox_inches='tight',pad_inches=0)

xs, dat = pool_data(all_avg_pks_cnt, timestamps, groups, pool_times, fill_missing)
plt.figure(figsize=(7,7))
plt.fill_between(range(1, len(xs)+1), [mean(e) - std(e) for e in dat],
                 [mean(e) + std(e) for e in dat], alpha=0.3)
plt.plot(range(1, len(xs)+1), [mean(e) for e in dat])
plt.xticks(ticks=range(1, len(dat) + 1), labels=xs, rotation=90)
plt.ylabel('Number of peaks')
plt.savefig(join(out_dir, "peaks_count.png"), dpi=300, bbox_inches='tight',pad_inches=0)

with open(join(out_dir, "peaks_count.txt"), 'w') as f:
    f.write("Time, AVG, STD\n")
    for i in range(len(xs)):
        f.write("{}, {}, {}\n".format(xs[i], mean(dat[i]), std(dat[i])))

plt.figure(figsize=(7,7))
dat = pool_data_avg_per_well(all_avg_pks_cnt, list(groups.keys())[0], timestamps, fill_missing)
for k in all_avg_pks_amps[list(groups.keys())[0]].keys():
    plt.errorbar(dat[k][0], [mean(e) for e in dat[k][1]], [std(e) for e in dat[k][1]])
plt.xticks(ticks=range(1, len(timestamps) + 1), labels=timestamps, rotation=90)
plt.ylabel('Number of peaks')
#plt.ylim(ymin = 0)
plt.savefig(join(out_dir, "peaks_count_wells.png"), dpi=300, bbox_inches='tight',pad_inches=0)

with open(join(out_dir, "peaks_count_wells.txt"), 'w') as f:
    f.write("Time, {}\n".format(", ".join(["Count_" + k for k in dat.keys()])))
    for i in range(len(xs)):
        vals = []
        for k in dat:
            if i < len(dat[k][1]):
                vals.append(str(dat[k][1][i][0]))
            else:
                vals.append("NaN")
        f.write("{}, {}\n".format(xs[i], ", ".join(vals)))

if time_groups:
    dats_t = pool_data_timecat(all_avg_pks_cnt, timestamps, list(groups.keys())[0], time_groups, fill_missing)
    plt.figure(figsize=(7,7))
    plt.bar([time_groups_names[k] for k in set(time_groups)], [mean(dats_t[k]) for k in set(time_groups)])
    plt.ylabel('Number of peaks')
    plt.savefig(join(out_dir, "peaks_count_bar.png"), dpi=300, bbox_inches='tight',pad_inches=0)

    with open(join(out_dir, "peaks_count_time.txt"), 'w') as f:
        f.write(", ".join([str(k) for k in set(time_groups)]) + "\n")
        for i in range(max([len(d) for d in dats_t.values()])):
            f.write(", ".join(["{:.3f}".format(dats_t[k][i]) for k in dats_t.keys()]) + "\n")

    print("Peak count:")
    for i in set(time_groups):
        print("{}: AVG ± SD = {:.2f} ± {:.2f}".format(time_groups_names[i], mean(dats_t[i]), std(dats_t[i])))
    for i in set(time_groups):
        for j in set(time_groups):
            if j > i:
                print("ANOVA {} vs {}: p = {:.6f}"
                      .format(time_groups_names[i], time_groups_names[j], f_oneway(dats_t[i], dats_t[j]).pvalue))


xs, dat = pool_data(all_traces_cnt, timestamps, groups, pool_times, fill_missing)
plt.figure(figsize=(7,7))
plt.fill_between(range(1, len(xs)+1), [mean(e) -std(e) for e in dat],
                 [mean(e) + std(e) for e in dat], alpha=0.3)
plt.plot(range(1, len(xs)+1), [mean(e) for e in dat])
plt.xticks(ticks=range(1, len(dat) + 1), labels=xs, rotation=90)
plt.ylabel('Number of traces')
plt.savefig(join(out_dir, "traces_count.png"), dpi=300, bbox_inches='tight',pad_inches=0)

with open(join(out_dir, "traces_count.txt"), 'w') as f:
    f.write("Time, AVG, STD\n")
    for i in range(len(xs)):
        f.write("{}, {}, {}\n".format(xs[i], mean(dat[i]), std(dat[i])))

plt.figure(figsize=(7,7))
dat = pool_data_avg_per_well(all_traces_cnt, list(groups.keys())[0], timestamps, fill_missing)
for k in dat.keys():
    plt.errorbar(dat[k][0], [mean(e) for e in dat[k][1]], [std(e) for e in dat[k][1]])
plt.xticks(ticks=range(1, len(timestamps) + 1), labels=timestamps, rotation=90)
plt.ylabel('Number of traces')
#plt.ylim(ymin = 0)
plt.savefig(join(out_dir, "traces_count_wells.png"), dpi=300, bbox_inches='tight',pad_inches=0)

with open(join(out_dir, "traces_count_wells.txt"), 'w') as f:
    f.write("Time, {}\n".format(", ".join(["Count_" + k for k in dat.keys()])))
    for i in range(len(xs)):
        vals = []
        for k in dat:
            if i < len(dat[k][1]):
                vals.append(str(dat[k][1][i][0]))
            else:
                vals.append("NaN")
        f.write("{}, {}\n".format(xs[i], ", ".join(vals)))

if time_groups:
    dats_t = pool_data_timecat(all_traces_cnt, timestamps, list(groups.keys())[0], time_groups, fill_missing)
    plt.figure(figsize=(7,7))
    plt.bar([time_groups_names[k] for k in set(time_groups)], [mean(dats_t[k]) for k in set(time_groups)])
    plt.ylabel('Number of traces')
    plt.savefig(join(out_dir, "traces_count_bar.png"), dpi=300, bbox_inches='tight',pad_inches=0)

    with open(join(out_dir, "traces_count_time.txt"), 'w') as f:
        f.write(", ".join([str(k) for k in set(time_groups)]) + "\n")
        for i in range(max([len(d) for d in dats_t.values()])):
            f.write(", ".join(["{:.3f}".format(dats_t[k][i]) for k in dats_t.keys()]) + "\n")

    print("traces count:")
    for i in set(time_groups):
        print("{}: AVG ± SD = {:.2f} ± {:.2f}".format(time_groups_names[i], mean(dats_t[i]), std(dats_t[i])))
    for i in set(time_groups):
        for j in set(time_groups):
            if j > i:
                print("ANOVA {} vs {}: p = {:.6f}"
                      .format(time_groups_names[i], time_groups_names[j], f_oneway(dats_t[i], dats_t[j]).pvalue))



# xs, dat = pool_data(all_avg_pks_freqs, timestamps, groups, pool_times)
# v_x = []
# v_y = []
# for i in range(len(dat)):
#     v_x.extend([i] * len(dat[i]))
#     v_y.extend(dat[i])
# plt.figure(figsize=(7,7))
# p = sns.violinplot(x=v_x, y=v_y, cut=0)
# plt.ylabel('Peak Frequency (Hz)')
# p.set_xticklabels(timestamps, rotation=90)
# plt.ylim([0, 0.4])
# plt.savefig(join(out_dir, "peaks_freq.png"), dpi=300, bbox_inches='tight',pad_inches=0)

xs, dat = pool_data(all_avg_pks_freqs, timestamps, groups, pool_times, fill_missing)
plt.figure(figsize=(7,7))
plt.fill_between(range(1, len(xs)+1), [mean(e) -std(e) for e in dat],
                 [mean(e) + std(e) for e in dat], alpha=0.3)
plt.plot(range(1, len(xs)+1), [mean(e) for e in dat])
plt.xticks(ticks=range(1, len(dat) + 1), labels=xs, rotation=90)
#plt.ylim([0, 0.8])
plt.ylabel('Peak Frequency (Hz)')
plt.savefig(join(out_dir, "peaks_freq.png"), dpi=300, bbox_inches='tight',pad_inches=0)

with open(join(out_dir, "peaks_freq.txt"), 'w') as f:
    f.write("Time, AVG, STD\n")
    for i in range(len(xs)):
        f.write("{}, {}, {}\n".format(xs[i], mean(dat[i]), std(dat[i])))


plt.figure(figsize=(7,7))
dat = pool_data_avg_per_well(all_avg_pks_freqs, list(groups.keys())[0], timestamps, fill_missing)
for k in dat.keys():
    plt.errorbar(dat[k][0], [mean(e) for e in dat[k][1]], [std(e) for e in dat[k][1]])
plt.xticks(ticks=range(1, len(timestamps) + 1), labels=timestamps, rotation=90)
plt.ylabel('Peak Frequency (Hz)')
plt.ylim(ymin = 0)
plt.savefig(join(out_dir, "peaks_freq_wells.png"), dpi=300, bbox_inches='tight',pad_inches=0)

with open(join(out_dir, "peaks_freq_wells.txt"), 'w') as f:
    f.write("Time, {}\n".format(", ".join(["AVG_" + k + ", STD_" + k for k in dat.keys()])))
    for i in range(len(xs)):
        vals = []
        for k in dat:
            if i < len(dat[k][1]):
                vals.append("{:.3f}, {:.3f}".format(mean(dat[k][1][i]), std(dat[k][1][i])))
            else:
                vals.append("NaN, NaN")
        f.write("{}, {}\n".format(xs[i], ", ".join(vals)))


if time_groups:
    dats_t = pool_data_timecat(all_avg_pks_freqs, timestamps, list(groups.keys())[0], time_groups, fill_missing)
    plt.figure(figsize=(7,7))
    plt.bar([time_groups_names[k] for k in set(time_groups)], [mean(dats_t[k]) for k in set(time_groups)])
    plt.ylabel('Peak Frequency (Hz)')
    plt.savefig(join(out_dir, "peaks_freq_bar.png"), dpi=300, bbox_inches='tight',pad_inches=0)

    with open(join(out_dir, "peaks_freq_time.txt"), 'w') as f:
        f.write(", ".join([str(k) for k in set(time_groups)]) + "\n")
        for i in range(max([len(d) for d in dats_t.values()])):
            f.write(", ".join(["{:.3f}".format(dats_t[k][i]) for k in dats_t.keys()]) + "\n")

    print("Peak frequency:")
    for i in set(time_groups):
        print("{}: AVG ± SD = {:.2f} ± {:.2f}".format(time_groups_names[i], mean(dats_t[i]), std(dats_t[i])))
    for i in set(time_groups):
        for j in set(time_groups):
            if j > i:
                print("ANOVA {} vs {}: p = {:.6f}"
                      .format(time_groups_names[i], time_groups_names[j], f_oneway(dats_t[i], dats_t[j]).pvalue))

#plt.figure(figsize=(7,7))
#p = plt.boxplot(dat)
#for i in range(len(dat)):
#    plt.plot(i + 1 + (rand(len(dat[i])) - 0.5) * 0.2, dat[i], 'x' + groups_col[gkeys[i % len(groups.keys())]])
#plt.xticks(ticks=range(1, len(dat) + 1), labels=xs, rotation=90)
#plt.ylabel('Peak Frequency (Hz)')
#plt.savefig(join(out_dir, "peaks_freq.png"), dpi=300, bbox_inches='tight',pad_inches=0)


# xs, dat = pool_data(all_avg_pks_amps, timestamps, groups, pool_times)
# v_x = []
# v_y = []
# for i in range(len(dat)):
#     v_x.extend([i] * len(dat[i]))
#     v_y.extend(dat[i])
# plt.figure(figsize=(7,7))
# p = sns.violinplot(x=v_x, y=v_y, cut=0)
# plt.ylabel('Peak amplitude (AU)')
# p.set_xticklabels(timestamps, rotation=90)
# plt.savefig(join(out_dir, "peaks_amp.png"), dpi=300, bbox_inches='tight',pad_inches=0)

# xs, dat = pool_data(all_avg_pks_amps, timestamps, groups, pool_times)
# plt.figure(figsize=(7,7))
# plt.fill_between(range(1, len(xs)+1), [mean(e) -std(e) for e in dat],
#                  [mean(e) + std(e) for e in dat], alpha=0.3)
# plt.plot(range(1, len(xs)+1), [mean(e) for e in dat])
# plt.xticks(ticks=range(1, len(dat) + 1), labels=xs, rotation=90)
# plt.ylabel('Peak amplitude (AU)')
#plt.savefig(join(out_dir, "peaks_freq.png"), dpi=300, bbox_inches='tight',pad_inches=0)

xs, dat = pool_data(all_avg_pks_amps, timestamps, groups, pool_times, fill_missing)
plt.figure(figsize=(7,7))
plt.fill_between(range(1, len(xs)+1), [mean(e) - std(e) for e in dat],
                 [mean(e) + std(e) for e in dat], alpha=0.3)
plt.plot(range(1, len(xs)+1), [mean(e) for e in dat])
plt.xticks(ticks=range(1, len(dat) + 1), labels=xs, rotation=90)
plt.ylabel('Peak amplitude (AU)')
plt.savefig(join(out_dir, "peaks_amplitude.png"), dpi=300, bbox_inches='tight',pad_inches=0)

with open(join(out_dir, "peaks_amplitude.txt"), 'w') as f:
    f.write("Time, AVG, STD\n")
    for i in range(len(xs)):
        f.write("{}, {}, {}\n".format(xs[i], mean(dat[i]), std(dat[i])))


plt.figure(figsize=(7,7))
dat = pool_data_avg_per_well(all_avg_pks_amps, list(groups.keys())[0], timestamps, fill_missing)
for k in dat.keys():
    plt.errorbar(dat[k][0], [mean(e) for e in dat[k][1]], [std(e) for e in dat[k][1]])
plt.xticks(ticks=range(1, len(timestamps) + 1), labels=timestamps, rotation=90)
plt.ylabel('Peak amplitude (AU)')
plt.ylim(ymin = 0)
plt.savefig(join(out_dir, "peaks_amplitude_wells.png"), dpi=300, bbox_inches='tight',pad_inches=0)

with open(join(out_dir, "peaks_amplitude_wells.txt"), 'w') as f:
    f.write("Time, {}\n".format(", ".join(["AVG_" + k + ", STD_" + k for k in dat.keys()])))
    for i in range(len(xs)):
        vals = []
        for k in dat:
            if i < len(dat[k][1]):
                vals.append("{:.3f}, {:.3f}".format(mean(dat[k][1][i]), std(dat[k][1][i])))
            else:
                vals.append("NaN, NaN")
        f.write("{}, {}\n".format(xs[i], ", ".join(vals)))

if time_groups:
    dats_t = pool_data_timecat(all_avg_pks_amps, timestamps, list(groups.keys())[0], time_groups, fill_missing)
    plt.figure(figsize=(7,7))
    plt.bar([time_groups_names[k] for k in set(time_groups)], [mean(dats_t[k]) for k in set(time_groups)])
    plt.ylabel('Peak amplitude (AU)')
    plt.savefig(join(out_dir, "peaks_amplitude_bar.png"), dpi=300, bbox_inches='tight',pad_inches=0)

    with open(join(out_dir, "peaks_amplitude_time.txt"), 'w') as f:
        f.write(", ".join([str(k) for k in set(time_groups)]) + "\n")
        for i in range(max([len(d) for d in dats_t.values()])):
            f.write(", ".join(["{:.3f}".format(dats_t[k][i]) for k in dats_t.keys()]) + "\n")

    print("Peak amplitude:")
    for i in set(time_groups):
        print("{}: AVG ± SD = {:.2f} ± {:.2f}".format(time_groups_names[i], mean(dats_t[i]), std(dats_t[i])))
    for i in set(time_groups):
        for j in set(time_groups):
            if j > i:
                print("ANOVA {} vs {}: p = {:.6f}"
                      .format(time_groups_names[i], time_groups_names[j], f_oneway(dats_t[i], dats_t[j]).pvalue))

#plt.figure(figsize=(7,7))
#p = plt.boxplot(dat)
#for i in range(len(dat)):
#    plt.plot(i + 1 + (rand(len(dat[i])) - 0.5) * 0.2, dat[i], 'x' + groups_col[gkeys[i % len(groups.keys())]])
#plt.xticks(ticks=range(1, len(dat) + 1), labels=xs, rotation=90)
#plt.ylabel('Peak amplitude (AU)')
#plt.savefig(join(out_dir, "peaks_amp.png"), dpi=300, bbox_inches='tight',pad_inches=0)


# xs, dat = pool_data(all_avg_pks_fwhm, timestamps, groups, pool_times)
# v_x = []
# v_y = []
# for i in range(len(dat)):
#     v_x.extend([i] * len(dat[i]))
#     v_y.extend(dat[i])
# plt.figure(figsize=(7,7))
# p = sns.violinplot(x=v_x, y=v_y, cut=0)
# p.set_xticklabels(timestamps, rotation=90)
# plt.ylabel('Peak width (s)')
# plt.savefig(join(out_dir, "peaks_width.png"), dpi=300, bbox_inches='tight',pad_inches=0)

xs, dat = pool_data(all_avg_pks_fwhm, timestamps, groups, pool_times, fill_missing)
plt.figure(figsize=(7,7))
plt.fill_between(range(1, len(xs)+1), [mean(e) -std(e) for e in dat],
                 [mean(e) + std(e) for e in dat], alpha=0.3)
plt.plot(range(1, len(xs)+1), [mean(e) for e in dat])
plt.xticks(ticks=range(1, len(dat) + 1), labels=xs, rotation=90)
plt.ylabel('Peak width (s)')
plt.savefig(join(out_dir, "peaks_width.png"), dpi=300, bbox_inches='tight',pad_inches=0)

with open(join(out_dir, "peaks_width.txt"), 'w') as f:
    f.write("Time, AVG, STD\n")
    for i in range(len(xs)):
        f.write("{}, {}, {}\n".format(xs[i], mean(dat[i]), std(dat[i])))

plt.figure(figsize=(7,7))
dat = pool_data_avg_per_well(all_avg_pks_fwhm, list(groups.keys())[0], timestamps, fill_missing)
for k in dat.keys():
    plt.errorbar(dat[k][0], [mean(e) for e in dat[k][1]], [std(e) for e in dat[k][1]])
plt.xticks(ticks=range(1, len(timestamps) + 1), labels=timestamps, rotation=90)
plt.ylabel('Peak width (s)')
plt.ylim(ymin = 0)
plt.savefig(join(out_dir, "peaks_width_wells.png"), dpi=300, bbox_inches='tight',pad_inches=0)


with open(join(out_dir, "peaks_width_wells.txt"), 'w') as f:
    f.write("Time, {}\n".format(", ".join(["AVG_" + k + ", STD_" + k for k in dat.keys()])))
    for i in range(len(xs)):
        vals = []
        for k in dat:
            if i < len(dat[k][1]):
                vals.append("{:.3f}, {:.3f}".format(mean(dat[k][1][i]), std(dat[k][1][i])))
            else:
                vals.append("NaN, NaN")
        f.write("{}, {}\n".format(xs[i], ", ".join(vals)))

if time_groups:
    dats_t = pool_data_timecat(all_avg_pks_fwhm, timestamps, list(groups.keys())[0], time_groups, fill_missing)
    plt.figure(figsize=(7,7))
    plt.bar([time_groups_names[k] for k in set(time_groups)], [mean(dats_t[k]) for k in set(time_groups)])
    plt.ylabel('AVG Peak width (s)')
    plt.savefig(join(out_dir, "peaks_width_bar.png"), dpi=300, bbox_inches='tight',pad_inches=0)

    with open(join(out_dir, "peaks_width_time.txt"), 'w') as f:
        f.write(", ".join([str(k) for k in set(time_groups)]) + "\n")
        for i in range(max([len(d) for d in dats_t.values()])):
            f.write(", ".join(["{:.3f}".format(dats_t[k][i]) for k in dats_t.keys()]) + "\n")

    print("Peak width:")
    for i in set(time_groups):
        print("{}: AVG ± SD = {:.2f} ± {:.2f}".format(time_groups_names[i], mean(dats_t[i]), std(dats_t[i])))
    for i in set(time_groups):
        for j in set(time_groups):
            if j > i:
                print("ANOVA {} vs {}: p = {:.6f}"
                      .format(time_groups_names[i], time_groups_names[j], f_oneway(dats_t[i], dats_t[j]).pvalue))


xs, dat = pool_data(all_max_ccor, timestamps, groups, pool_times, fill_missing)
plt.figure(figsize=(7,7))
plt.fill_between(range(1, len(xs)+1), [mean(e) - std(e) for e in dat],
                 [mean(e) + std(e) for e in dat], alpha=0.3)
plt.plot(range(1, len(xs)+1), [mean(e) for e in dat])
plt.xticks(ticks=range(1, len(dat) + 1), labels=xs, rotation=90)
plt.ylabel('AVG max correlation')
#plt.ylim([0.6, 1])
plt.savefig(join(out_dir, "corr.png"), dpi=300, bbox_inches='tight',pad_inches=0)

with open(join(out_dir, "corr.txt"), 'w') as f:
    f.write("Time, AVG, STD\n")
    for i in range(len(xs)):
        f.write("{}, {:.3f}, {:.3f}\n".format(xs[i], mean(dat[i]), std(dat[i])))

plt.figure(figsize=(7,7))
dat = pool_data_avg_per_well(all_max_ccor, list(groups.keys())[0], timestamps, fill_missing)
for k in dat.keys():
    plt.errorbar(dat[k][0], [mean(e) for e in dat[k][1]], [std(e) for e in dat[k][1]])
plt.xticks(ticks=range(1, len(timestamps) + 1), labels=timestamps, rotation=90)
plt.ylabel('AVG max correlation')
plt.ylim(ymax = 1)
plt.savefig(join(out_dir, "corr_wells.png"), dpi=300, bbox_inches='tight',pad_inches=0)

with open(join(out_dir, "corr_wells.txt"), 'w') as f:
    f.write("Time, {}\n".format(", ".join(["AVG_" + k + ", STD_" + k for k in dat.keys()])))
    for i in range(len(xs)):
        vals = []
        for k in dat:
            if i < len(dat[k][1]):
                vals.append("{:.3f}, {:.3f}".format(mean(dat[k][1][i]), std(dat[k][1][i])))
            else:
                vals.append("NaN, NaN")
        f.write("{}, {}\n".format(xs[i], ", ".join(vals)))

if time_groups:
    dats_t = pool_data_timecat(all_max_ccor, timestamps, list(groups.keys())[0], time_groups, fill_missing)
    plt.figure(figsize=(7,7))
    plt.bar([time_groups_names[k] for k in set(time_groups)], [mean(dats_t[k]) for k in set(time_groups)])
    plt.ylabel('AVG max correlation')
    plt.savefig(join(out_dir, "corr_bar.png"), dpi=300, bbox_inches='tight',pad_inches=0)

    with open(join(out_dir, "corr_time.txt"), 'w') as f:
        f.write(", ".join([str(k) for k in set(time_groups)]) + "\n")
        for i in range(max([len(d) for d in dats_t.values()])):
            f.write(", ".join(["{:.3f}".format(dats_t[k][i]) for k in dats_t.keys()]) + "\n")

    print("Peak correlation:")
    for i in set(time_groups):
        print("{}: AVG ± SD = {:.2f} ± {:.2f}".format(time_groups_names[i], mean(dats_t[i]), std(dats_t[i])))
    for i in set(time_groups):
        for j in set(time_groups):
            if j > i:
                print("ANOVA {} vs {}: p = {:.6f}"
                      .format(time_groups_names[i], time_groups_names[j], f_oneway(dats_t[i], dats_t[j]).pvalue))
