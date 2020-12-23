#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  1 19:19:50 2020

@author: pierre
"""

from os.path import join, isfile
import pickle

from numpy import array, mean, max, std, histogram, correlate, argmax, arange, sqrt, cumsum
from numpy.polynomial.polynomial import polyfit

from matplotlib import pyplot as plt
from scipy.signal import find_peaks


from scipy.stats import pearsonr

from math import isnan

cache_dir = "/mnt/data2/mosab_incucyte/processed/"

DT = 0.33
pxsize = 0.2646 #mm

well_id = "B3_1"
timestamps = ["01d19h48m", "02d10h48m", "02d20h21m"]

colors = {timestamps[0]:"#00FF49FF", timestamps[1]: "#FF4900FF", timestamps[2]: "#4900FFFF"}
colors_alpha = {timestamps[0]:"#00FF494B", timestamps[1]: "#FF49004B", timestamps[2]: "#4900FF4B"}


traces = {}
objs_cents = {}
for timestamp in timestamps:
    fname = "/mnt/data2/mosab_incucyte/processed/traces_" + well_id + "_" + timestamp + ".csv"
    traces[timestamp] = []
    with open(fname, 'r') as f:
        for line in f:
            line = line.rstrip("\n").split(" ")
            traces[timestamp].append([float(e) for e in line[1].split(",")])

    fname = "/mnt/data2/mosab_incucyte/processed/objects_" + well_id + "_" + timestamp + ".csv"
    objs_cents[timestamp] = []
    with open(fname, 'r') as f:
        for line in f:
            line = line.strip("\n").split(" ")
            tmp = array([[float(e.split(",")[0]), float(e.split(",")[1])] for e in line[1:]])
            objs_cents[timestamp].append(mean(tmp, axis=0))

avg_sig = {timestamp: [mean([e[i] for e in traces[timestamp]]) for i in range(len(traces[timestamp][0]))]
            for timestamp in timestamps}
avg_sig_pks =  {ts:[e for e in find_peaks(avg_sig[ts])[0] if avg_sig[ts][e] > 20] for ts in timestamps}

plt.figure()
plt.bar(timestamps, [len(traces[ts]) for ts in timestamps])
plt.ylabel('Number of detected neurons')
plt.savefig("/tmp/num_neurons" + well_id+ ".png", dpi=300, bbox_inches='tight')

print("Num cells")
for ts in timestamps:
    print("{}: {}".format(ts, len(traces[ts])))

for ts in timestamps:
    plt.figure()
    for trace in traces[ts]:
        plt.plot(array(range(len(trace))) * DT, trace, linewidth=0.5)
    plt.plot(array(range(len(avg_sig[ts]))) * DT, avg_sig[ts], 'k', linewidth=2)
    #plt.plot([e * DT for e in avg_sig_pks[ts]], [avg_sig[ts][e] for e in avg_sig_pks[ts]], 'xr')
    plt.ylim([0, 115])
    plt.ylabel("Intensity (AU)")
    plt.xlabel("Time (s)")
    plt.savefig("/tmp/traces_" + well_id + "_" + ts + ".png", dpi=300, bbox_inches='tight')

print("Trace intensities")
for ts in timestamps:
    print("{}: {:.2f} ± {:.2f}".format(ts, mean([max(trace) for trace in traces[ts]]), std([max(trace) for trace in traces[ts]])))

plt.figure()
o, b = histogram([max(trace) for trace in traces[timestamps[0]]], bins=30)
plt.bar(b[:-1], o / sum(o), color="#00FF494B", align="edge", width=b[1]- b[0])
o, b = histogram([max(trace) for trace in traces[timestamps[1]]], bins=30)
plt.bar(b[:-1], o / sum(o), color="#FF49004B", align="edge", width=b[1]- b[0])
o, b = histogram([max(trace) for trace in traces[timestamps[2]]], bins=30)
plt.bar(b[:-1], o / sum(o), color="#4900FF4B", align="edge", width=b[1]- b[0])
plt.xlabel('Max cell intensity (AU)')
plt.ylabel('Frequency')
plt.legend(timestamps)
plt.savefig("/tmp/intensity_hists_" + well_id+ ".png", dpi=300, bbox_inches='tight')


print("Mean Trace intensities")
for ts in timestamps:
    print("{}: {:.2f} ± {:.2f}".format(ts, mean([mean(trace) for trace in traces[ts]]), std([mean(trace) for trace in traces[ts]])))

plt.figure()
for ts in timestamps:
    o, b = histogram([mean(trace) for trace in traces[ts]], bins=30)
    plt.bar(b[:-1], o / sum(o), color=colors_alpha[ts], align="edge", width=b[1]- b[0])
plt.xlabel('Mean cell intensity (AU)')
plt.ylabel('Frequency')
plt.legend(timestamps)
plt.savefig("/tmp/mean_intensity_hists_" + well_id+ ".png", dpi=300, bbox_inches='tight')


tmp = {ts:[] for ts in timestamps}
[[tmp[ts].extend([e[f] for f in avg_sig_pks[ts]]) for e in traces[ts]] for ts in timestamps]

print("Intensities at avg signal peaks")
for ts in timestamps:
    print("{}: {:.2f} ± {:.2f}".format(ts, mean(tmp[ts]), std(tmp[ts])))

plt.figure()
for ts in timestamps:
    o, b = histogram(tmp[ts], bins=30)
    plt.bar(b[:-1], o / sum(o), color=colors_alpha[ts], align="edge", width=b[1]- b[0])
plt.xlabel('Mean intensities at avg signal peaks (AU)')
plt.ylabel('Frequency')
plt.legend(timestamps)
plt.savefig("/tmp/mean_intensity_at_avg_peak_hists_" + well_id+ ".png", dpi=300, bbox_inches='tight')




plt.figure()
o, b = histogram([max(trace) for trace in traces[timestamps[0]]], arange(0, 120, 2))
plt.plot(b[:-1], cumsum(o), color="#00FF494B")
o, b = histogram([max(trace) for trace in traces[timestamps[1]]], arange(0, 120, 2))
plt.plot(b[:-1], cumsum(o), color="#FF49004B")
o, b = histogram([max(trace) for trace in traces[timestamps[2]]], arange(0, 120, 2))
plt.plot(b[:-1], cumsum(o), color="#4900FF4B")
plt.xlabel('Max cell intensity (AU)')
plt.ylabel('Number of cells')
plt.legend(timestamps)
plt.savefig("/tmp/intensity_cum_" + well_id+ ".png", dpi=300, bbox_inches='tight')


plt.figure()
for ts in timestamps:
    pks = find_peaks(avg_sig[ts])[0]
    plt.plot(array(range(0, len(avg_sig[ts]) - pks[0])) * DT, avg_sig[ts][pks[0]:], color=colors[ts])
plt.legend(timestamps)
plt.ylabel("Intensity (AU)")
plt.xlabel("Time to first peak (s)")
plt.ylim([0, 55])
#plt.xlim([0, 175])
plt.savefig("/tmp/intensity_averages_" + well_id+ ".png", dpi=300, bbox_inches='tight')

pks_delta = {}
pks_amp = {}
for timestamp in timestamps:
    pks_delta[timestamp] = []
    pks_amp[timestamp] = []
    for trace in traces[timestamp]:
        cur_pks = find_peaks(trace)[0]
        Mtr = mean(trace)
        SDtr = std(trace)

        cur_pks = [pk for pk in cur_pks if trace[pk] > 2*Mtr]

        pks_delta[timestamp].append([(cur_pks[i+1] - cur_pks[i]) * DT for i in range(0, len(cur_pks) - 1) if (cur_pks[i+1] - cur_pks[i]) * DT < 50])
        pks_amp[timestamp].append([trace[pk] for pk in cur_pks])
plt.savefig("/tmp/yoyo" + well_id+ ".png", dpi=300, bbox_inches='tight')




tmp = {ts:[] for ts in timestamps}
[[tmp[ts].extend(e) for e in pks_delta[ts]] for ts in timestamps]

fig = plt.figure()
o, b = histogram(tmp[timestamps[0]], bins=array(range(0, 120, 4)) * DT)
plt.bar(b[:-1], o / sum(o), color="#00FF494B", align="edge", width=b[1]- b[0])
o, b = histogram(tmp[timestamps[1]], bins=array(range(0, 120, 4)) * DT)
plt.bar(b[:-1], o / sum(o), color="#FF49004B", align="edge", width=b[1]- b[0])
o, b = histogram(tmp[timestamps[2]], bins=array(range(0, 120, 4)) * DT)
plt.bar(b[:-1], o / sum(o), color="#4900FF4B", align="edge", width=b[1]- b[0])
plt.xlabel('Time between peaks (s)')
plt.ylabel('Frequency')
plt.legend(timestamps)
plt.savefig("/tmp/peak_period_" + well_id+ ".png", dpi=300, bbox_inches='tight')

print("Peak period")
for ts in timestamps:
    print("{}: {:.2f} ± {:.2f}".format(ts, mean(tmp[ts]), std(tmp[ts])))

print("Peak frequency")
for ts in timestamps:
    print("{}: {:.2f} ± {:.2f}".format(ts, mean([1 / e for e in tmp[ts]]), std([1 / e for e in tmp[ts]])))

tmp = {ts:[] for ts in timestamps}
[[tmp[ts].extend(e) for e in pks_amp[ts]] for ts in timestamps]

plt.figure()
o, b = histogram(tmp[timestamps[0]], bins=40)
plt.bar(b[:-1], o / sum(o), color="#00FF494B", align="edge", width=b[1]- b[0])
o, b = histogram(tmp[timestamps[1]], bins=40)
plt.bar(b[:-1], o / sum(o), color="#FF49004B", align="edge", width=b[1]- b[0])
o, b = histogram(tmp[timestamps[2]], bins=40)
plt.bar(b[:-1], o / sum(o), color="#4900FF4B", align="edge", width=b[1]- b[0])
plt.xlabel('Peak intensity (AU)')
plt.ylabel("Frequency")
plt.legend(timestamps)
plt.savefig("/tmp/peak_intens_" + well_id+ ".png", dpi=300, bbox_inches='tight')

print("Peak intensity")
for ts in timestamps:
    print("{}: {:.2f} ± {:.2f}".format(ts, mean(tmp[ts]), std(tmp[ts])))


corr_cache_f = join(cache_dir, "correlations_" + well_id + "_" + timestamp + ".pkl")
if isfile(corr_cache_f):
    with open(corr_cache_f, "rb") as f:
        dat = pickle.load(f)
        #rs = dat["rs"]
        ccs_max = dat["ccs_max"]
        ccs_0 = dat["ccs_0"]
        ccs_shift = dat["ccs_shift"]
        avg_ccs_max = dat["avg_ccs_max"]
else:
    #rs = {ts:[] for ts in timestamps}
    #ccs = {ts:[] for ts in timestamps}
    ccs_0 = {ts:[] for ts in timestamps}
    ccs_max = {ts:[] for ts in timestamps}
    ccs_shift = {ts:[] for ts in timestamps}
    avg_ccs_max = {ts:[] for ts in timestamps}
    for ts in timestamps:
        for i in range(len(traces[ts])): #range(1):
            #rs[ts].append([])
            #ccs[ts].append([])
            ccs_0[ts].append([])
            ccs_max[ts].append([])
            ccs_shift[ts].append([])
            tmp = []
            for j in range(len(traces[ts])):
                if i == j:
                    #rs[ts][-1].append([])
                    continue
                #rs[ts][-1].append(pearsonr(traces[ts][i], traces[ts][j]))
                ccov = correlate(traces[ts][i] - mean(traces[ts][i]), traces[ts][j] - mean(traces[ts][j]), mode='full')
                ccor = ccov / (len(traces[ts][i]) * std(traces[ts][i]) * std(traces[ts][j]))
                #ccs[ts][-1].append(ccor)

                tmp.append(max(ccor))
                if j > i:
                    ccs_0[ts][-1].append(ccor[int((len(ccor) - 1) / 2)])
                    ccs_max[ts][-1].append(max(ccor))
                    ccs_shift[ts][-1].append(argmax(ccor))
            avg_ccs_max[ts].append(mean([e for e in tmp if not isnan(e)]))

    with open(corr_cache_f, 'wb') as f:
        pickle.dump({"timestamps": timestamps, "ccs_max": ccs_max,
                     "ccs_shift": ccs_shift, "avg_ccs_max": avg_ccs_max,
                     "ccs_0": ccs_0}, f) #"rs": rs,


#idx_max = [0]
#plt.figure()
#for i in range(1, len(traces[ts])):
#    plt.plot(range(len(ccs[ts][0][i])), ccs[ts][0][i])
#    idx_max.append(argmax(ccs[ts][0][i]) - (len(ccs[ts][0][i]) - 1) / 2)

#plt.figure()
#plt.plot(range(len(traces[ts][0])), traces[ts][0], 'k')
#for i in range(1, len(traces[ts])):
#    plt.plot(arange(len(traces[ts][j])) + idx_max[i], traces[ts][j])
#plt.savefig("/tmp/test" + well_id+ ".png", dpi=300)

#plt.figure()
#plt.plot(range(len(ccs[ts][0][1])), ccs[ts][0][1])


# tst = correlate(traces[ts][0] - mean(traces[ts][0]), traces[ts][0] - mean(traces[ts][0]), mode='full')
# tst = tst / (len(traces[ts][0]) * std(traces[ts][0]) * std(traces[ts][0]))

# plt.figure()
# plt.plot(range(len(tst)), tst)

#plt.figure()
#plt.hist(idx_max, bins=range(-230, 50), density=True)
#plt.savefig('/tmp/a.png', dpi=300)


tmp = {ts:[] for ts in timestamps}
for ts in timestamps:
    [tmp[ts].extend([f for e in ccs_max[ts] for f in e if not isnan(f)])]


print("Mean max correlations:")
for ts in timestamps:
    print("{}: {:.2f} ± {:.2f}".format(ts, mean(tmp[ts]), std(tmp[ts])))

plt.figure()
o, b = histogram(tmp[timestamps[0]], bins=arange(0, 1, 0.01))
plt.bar(b[:-1], o / sum(o), color="#00FF494B", align="edge", width=b[1]- b[0])
o, b = histogram(tmp[timestamps[1]], bins=arange(0, 1, 0.01))
plt.bar(b[:-1], o / sum(o), color="#FF49004B", align="edge", width=b[1]- b[0])
o, b = histogram(tmp[timestamps[2]], bins=arange(0, 1, 0.01))
plt.bar(b[:-1], o / sum(o), color="#4900FF4B", align="edge", width=b[1]- b[0])
plt.xlabel('Correlation at optimal τ')
plt.ylabel('Frequency')
plt.legend(timestamps)
plt.savefig("/tmp/cross_corr_max_hist_" + well_id+ ".png", dpi=300, bbox_inches='tight')

tmp = {ts:[] for ts in timestamps}
for ts in timestamps:
    [tmp[ts].extend([f for e in ccs_0[ts] for f in e if not isnan(f)])]


print("Mean correlations at τ=0:")
for ts in timestamps:
    print("{}: {:.2f} ± {:.2f}".format(ts, mean(tmp[ts]), std(tmp[ts])))


plt.figure()
o, b = histogram(tmp[timestamps[0]], bins=arange(0, 1, 0.01))
plt.bar(b[:-1], o / sum(o), color="#00FF494B", align="edge", width=b[1]- b[0])
o, b = histogram(tmp[timestamps[1]], bins=arange(0, 1, 0.01))
plt.bar(b[:-1], o / sum(o), color="#FF49004B", align="edge", width=b[1]- b[0])
o, b = histogram(tmp[timestamps[2]], bins=arange(0, 1, 0.01))
plt.bar(b[:-1], o / sum(o), color="#4900FF4B", align="edge", width=b[1]- b[0])
plt.xlabel('Correlation at τ=0')
plt.ylabel('Frequency')
plt.legend(timestamps)
plt.savefig("/tmp/cross_corr_0_hist_" + well_id+ ".png", dpi=300, bbox_inches='tight')


plt.figure()
for ts in timestamps:
    o, b = histogram(tmp[ts], bins=arange(0, 1, 0.01))
    plt.bar(b[:-1], o, color=colors_alpha[ts], align="edge", width=b[1]- b[0])
plt.xlabel('Correlation at τ=0')
plt.ylabel('Count')
plt.legend(timestamps)
plt.savefig("/tmp/cross_corr_0_hist_cnt_" + well_id+ ".png", dpi=300, bbox_inches='tight')

print("Mean max avg correlations:")
for ts in timestamps:
    print("{}: {:.2f} ± {:.2f}".format(ts, mean(avg_ccs_max[ts]), std(avg_ccs_max[ts])))

plt.figure()
o, b = histogram(avg_ccs_max[timestamps[0]], bins=arange(0, 1, 0.01))
plt.bar(b[:-1], o / sum(o), color="#00FF494B", align="edge", width=b[1]- b[0])
o, b = histogram(avg_ccs_max[timestamps[1]], bins=arange(0, 1, 0.01))
plt.bar(b[:-1], o / sum(o), color="#FF49004B", align="edge", width=b[1]- b[0])
o, b = histogram(avg_ccs_max[timestamps[2]], bins=arange(0, 1, 0.01))
plt.bar(b[:-1], o / sum(o), color="#4900FF4B", align="edge", width=b[1]- b[0])
plt.xlabel('Average trace correlation')
plt.ylabel('Frequency')
plt.legend(timestamps)
plt.savefig("/tmp/cross_corr_avgmax_hist_" + well_id+ ".png", dpi=300, bbox_inches='tight')


tmp = {ts:[] for ts in timestamps}
for ts in timestamps:
    [tmp[ts].extend(e) for e in ccs_shift[ts]]

plt.figure()
o, b = histogram([e - (len(traces[timestamps[0]][0]) - 1) for e in tmp[timestamps[0]]], bins=arange(-20, 20, 1))
plt.bar(b[:-1], o / sum(o), color="#00FF494B", align="edge", width=b[1]- b[0])
o, b = histogram([e - (len(traces[timestamps[1]][0]) - 1) for e in tmp[timestamps[1]]], bins=arange(-20, 20, 1))
plt.bar(b[:-1], o / sum(o), color="#FF49004B", align="edge", width=b[1]- b[0])
o, b = histogram([e - (len(traces[timestamps[2]][0]) - 1) for e in tmp[timestamps[2]]], bins=arange(-20, 20, 1))
plt.bar(b[:-1], o / sum(o), color="#4900FF4B", align="edge", width=b[1]- b[0])
plt.xlabel('Cross-correlation phase shift')
plt.ylabel('Frequency')
plt.legend(timestamps)
plt.savefig("/tmp/cross_corr_shift_hist_" + well_id+ ".png", dpi=300, bbox_inches='tight')


print("Individual correlations")
idxs_ts = {timestamps[0]: 966, timestamps[1]: 968, timestamps[2]: 966}
dist_corrs = {}
for ts in timestamps:
    dist_corrs[ts] = []
    for i in range(len(avg_ccs_max[ts])):
        idxs = list(range(0, i)) + list(range(i+1, len(avg_ccs_max[ts])))
        xs = [sqrt(sum((objs_cents[ts][k] - objs_cents[ts][i])**2)) * pxsize for k in idxs]
        #ys = ccs_max[ts][i] + [ccs_max[ts][j][i-j-1] for j in range(i)]
        ys = [ccs_0[ts][j][i-j-1] for j in range(i)] + ccs_0[ts][i]
        dist_corrs[ts].append(pearsonr(xs, ys)[0])

        if i == idxs_ts[ts]:
            b, m = polyfit(xs, ys, 1)
            print("{},{}: r={}".format(ts, i, m))
            plt.figure()
            plt.plot(xs, ys, 'x', color=colors[ts])
            plt.plot(xs, b + m * array(xs), 'k')
            plt.ylabel('Max cross-correlation')
            plt.xlabel('Distance (mm)')
            plt.ylim([0, 1])
            plt.xlim([0, 300])
            plt.savefig("/tmp/cross_corr_dist_ex_" + ts + "_" + well_id + ".png", dpi=300, bbox_inches='tight')

print("Correlation distance")
for ts in timestamps:
    print("{}: {:.2f} ± {:.2f}".format(ts, mean(dist_corrs[ts]), std(dist_corrs[ts])))


plt.figure()
o, b = histogram(dist_corrs[timestamps[0]], bins=arange(-1, 1, 0.05))
plt.bar(b[:-1], o / sum(o), color="#00FF494B", align="edge", width=b[1]- b[0])
o, b = histogram(dist_corrs[timestamps[1]], bins=arange(-1, 1, 0.05))
plt.bar(b[:-1], o / sum(o), color="#FF49004B", align="edge", width=b[1]- b[0])
o, b = histogram(dist_corrs[timestamps[2]], bins=arange(-1, 1, 0.05))
plt.bar(b[:-1], o / sum(o), color="#4900FF4B", align="edge", width=b[1]- b[0])
plt.xlabel('Trace correlation vs distance')
plt.ylabel('Frequency')
plt.legend(timestamps)
plt.savefig("/tmp/cross_corr_dist_hist_" + well_id+ ".png", dpi=300, bbox_inches='tight')

corr_to_avg_0 = {ts:[] for ts in timestamps}
for ts in timestamps:
    norm_avg_sig = avg_sig[ts] - mean(avg_sig[ts])
    for i in range(len(traces[ts])):
        ccov = correlate(traces[ts][i] - mean(traces[ts][i]), norm_avg_sig, mode='full')
        ccor = ccov / (len(traces[ts][i]) * std(traces[ts][i]) * std(avg_sig[ts]))
        corr_to_avg_0[ts].append(ccor[int((len(ccor) - 1) / 2)])

print("Correlation with avg signal")
for ts in timestamps:
    print("{}: {:.2f} ± {:.2f}".format(ts, mean(corr_to_avg_0[ts]), std(corr_to_avg_0[ts])))

plt.figure()
for ts in timestamps:
    o, b = histogram(corr_to_avg_0[ts], bins=arange(0, 1, 0.01))
    plt.bar(b[:-1], o / sum(o), color=colors_alpha[ts], align="edge", width=b[1]- b[0])
plt.xlabel('Correlation to avg signal at τ=0')
plt.ylabel('Frequency')
plt.legend(timestamps)
plt.savefig("/tmp/cross_corr_0_avgsig_" + well_id+ ".png", dpi=300, bbox_inches='tight')
