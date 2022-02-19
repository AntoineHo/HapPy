#!/usr/bin/python
# -*- coding: utf-8 -*-

# General
import os, sys, math

# Happy
try :
    from happy.utils import *
except :
    from utils import *

try :
    from happy.plot import *
except :
    from plot import *

# Stats
from scipy.signal import savgol_filter  # SmoothingAUC
from scipy.signal import find_peaks  # Finding peaks
from scipy.signal import peak_widths

def auto_find_limits(peaks, thresholds) :
    """Given 2 peaks positions and thresholds, return a dict like
    {"low": int, "dip": int, "high": int} with limits to use"""

    limits = {}
    diploid_peak = min(peaks)
    haploid_peak = max(peaks)

    for t, lim in thresholds.items() :
        if diploid_peak in t and haploid_peak in t :
            limits["dip"] = lim
        elif diploid_peak in t and haploid_peak not in t :
            limits["low"] = lim
        elif haploid_peak in t and diploid_peak not in t :
            limits["high"] = lim
        else :
            continue

    if any(x not in limits for x in ["low", "dip", "high"]) :
        raise Exception("Could not automatically find limits...")

    return limits

def get_threshold_between_peaks(smoothed, peaks, valleys) :
    """For each pairs of consecutive peaks, find the most "middle" valley to set the threshold"""
    # For all pairs of consecutive : find the valley in between
    valleys = valleys
    tresholds = {}
    first_peak_width = peak_widths(smoothed, [peaks[0]])[0][0]  # Get peak widths
    first_peak_boundary = int(peaks[0] - first_peak_width)
    first_peak_boundary = 0 if first_peak_boundary < 0 else first_peak_boundary

    last_peak_width = peak_widths(smoothed, [peaks[-1]])[0][0]  # Get peak widths
    last_peak_boundary = int(peaks[-1] + last_peak_width)
    last_peak_boundary = len(smoothed) if last_peak_boundary > len(smoothed) else last_peak_boundary

    tresholds[(0,peaks[0])] = first_peak_boundary

    for p1, p2 in zip(peaks, peaks[1:]) :
        valid_thresholds = []
        for val in valleys :
            if p1 < val < p2 :
                valid_thresholds.append(val)
            else :
                continue

        if len(valid_thresholds) > 1 :

            most_middle_threshold, diff_size = None, None
            for v in valid_thresholds :
                x = v - p1
                y = p2 - v
                diff = abs(x - y)
                if diff_size is None or diff < diff_size :
                    most_middle_threshold = v
                    diff_size = diff
                else :
                    continue

            tresholds[(p1,p2)] = most_middle_threshold
        else :
            tresholds[(p1,p2)] = valid_thresholds[0]

    # last peak
    tresholds[(peaks[-1], "inf")] = last_peak_boundary

    return tresholds

def determine_peaks_and_limits(
    data, smoothed, limits,
    peak_prom, peak_height,
    valley_prom, valley_height,
    debug, smooth_window_size, outfile,
    skip_smooth,
):
    """Use a smoothed frequency histogram to find peaks and valleys
    Then use the determined peaks and valleys as limits for AUC computation
    """
    mm = max(smoothed)
    peaks, props = find_peaks(smoothed, height=peak_height, prominence=peak_prom) # maxima (peaks positions)
    rpeaks, rprops = find_peaks([-i+mm for i in smoothed], height=valley_height, prominence=valley_prom) # minima (peaks limits)

    if len(peaks) > 3 :
        print("WARNING: More than 3 peaks detected.\nPossible erroneous detection:\n\t-Restart setting the -ll parameter.\n\t-check histogram and modify peak height and prominence arguments accordingly.\n\t-Contaminant peaks may also break detection, remove them with tools such as blobtools or by hard-filtering low coverage contigs.")
        print("NOTE: Assuming the last 2 peaks are diploid and haploid...")

    if debug :
        debug_plot_peak_errors(data, smoothed, peaks, limits.values(), rpeaks, smooth_window_size, outfile, skip_smooth)

    if len(peaks) > 0 :
        print("Peaks found: " + "x, ".join(str(p) for p in peaks) + "x")
    else :
        raise Exception("No peaks found! Try changing the input parameters or setting thresholds manually!")
    if len(rpeaks) > 0 :
        print("Valleys found: " + "x, ".join(str(p) for p in rpeaks) + "x")
    else :
        print("No valleys found!")

    valleys = [0] + list(rpeaks) + [len(smoothed)]
    thresholds = get_threshold_between_peaks(smoothed, peaks, valleys)

    relevant_peaks = peaks[-3:]
    #valleys = rpeaks[-3:]
    print("Relevant peaks: " + "x, ".join(str(p) for p in relevant_peaks) + "x")
    print("Thresholds:\n\t- " + "\t- ".join("{}: {}x\n".format(k,p) for k,p in thresholds.items()))

    return relevant_peaks, thresholds

def estimate_haploidy(
    infile, size, outfile, # required arguments
    limit_low, limit_dip_hap, limit_high, # manual thresholds (optional)
    peak_prom, peak_height, # peak detection
    valley_prom, valley_height, # valley detection
    window, skip_smooth=False,
    plot=False, debug=False,
):
    """Finds peaks and modality, then computes scores of haploidy"""

    # Get histogram file and size estimation
    HIST = check_files(infile)
    SIZE = size_from_string(size)

    print("# Hap.py estimate")
    print("Coverage histogram:\t{0}\nOutput file:\t{1}\n".format(HIST, outfile))
    print(
        "===============================================================================\n"
    )

    # Read coverage histogram
    log("Reading histogram...")
    freqs = [int(line.split()[1]) for line in open(HIST, "r")][1:-1]

    # Apply Savitzky-Golay filter to smooth coverage histogram
    log("Analysing curve...")
    if skip_smooth :
        smoothed = freqs
    else :
        if len(freqs) <= window :
            print("WARNING: SKIPPED smoothing because max coverage < window size setting.\nTo avoid this warning set a lower '-w' or use '--no-smooth'.")
            smoothed = freqs
        else :
            smoothed = [s for s in savgol_filter(freqs, window, 3) if s >= 0]

    limits = {
        "low":  limit_low,
        "dip":  limit_dip_hap,
        "high": limit_high,
    }

    auto_detect = False
    if all(x is None for x in limits.values()) :
        log("Detecting peaks and thresholds automatically...")
        peaks, thresholds = determine_peaks_and_limits(
            freqs, smoothed, limits,
            peak_prom, peak_height, # peak detection
            valley_prom, valley_height, # valley detection
            debug, window, outfile, skip_smooth
        )
        auto_detect = True
    elif any(x is None for x in limits.values()) :
        raise Exception("You must set all values for thresholds -ll, -ld and -lh...")
    else :
        log("Using input thresholds to compute AUC...")
        #peaks, props = find_peaks(smoothed, height=peak_height, prominence=peak_prom) # maxima (peaks positions)
        peaks, props = find_peaks(smoothed[limits["low"] : limits["high"]], height=peak_height, prominence=peak_prom) # maxima (peaks positions)
        peaks = [p+limits["low"] for p in peaks]
        auto_detect = False

    heights = [smoothed[i] for i in peaks]

    if debug:
        debug_smooth_histogram(smoothed, peaks, heights, outfile, skip_smooth=skip_smooth)

    params, cov = None, None
    peak_ratios = None

    if not auto_detect :
        if len(peaks) > 2 :
            print("WARNING: detected more than 2 peaks in given limits: {}x to {}x".format(limits["low"], limits["high"]))
            print("WARNING: haplotigs peak ratio not computed...")
            haplotigs_peak_ratio = None
        elif len(peaks) < 2 :
            print("WARNING: detected less than 2 peaks in given limits: {}x to {}x".format(limits["low"], limits["high"]))
            print("WARNING: haplotigs peak ratio not computed...")
            haplotigs_peak_ratio = None
        else :
            # haplotigs_peak_ratio = diploid peak / haploid peak : if > 1 then more diploid if < 1 then more haploid
            haplotigs_peak_ratio = round(freqs[min(peaks)]/freqs[max(peaks)], 3)

        AUC_low = sum(freqs[: limits["low"]])
        AUC_diplo = sum(freqs[limits["low"] : limits["dip"]])
        AUC_haplo = sum(freqs[limits["dip"] : limits["high"]])
        AUC_high = sum(freqs[limits["high"]:])
    else :
        if len(peaks) == 0 : # Cannot compute
            raise Exception("No peak found! Try with --debug to check distribution.")
        elif len(peaks) == 1 : # If only one peak
            #log("Found 1 peak at: {}x".format(peaks[0]))
            log("WARNING: Only 1 peak found, either set thresholds manually with -ll, -ld and -lh or adapt peak detection parameters...")
            log("Finished!")
            sys.exit(0)

        elif len(peaks) == 2:  # If only 2 peaks
            #log("Assuming peaks are not low...: {0}x and {1}x".format(peaks[0], peaks[1]))
            limits = auto_find_limits(peaks, thresholds)
            AUC_low = sum(freqs[: limits["low"]])
            AUC_diplo = sum(freqs[limits["low"] : limits["dip"]])
            AUC_haplo = sum(freqs[limits["dip"] : limits["high"]])
            AUC_high = sum(freqs[limits["high"]:])
            haplotigs_peak_ratio = round(freqs[min(peaks)]/freqs[max(peaks)], 3)

        elif len(peaks) == 3:  # In case 3 peaks : 1 in each category
            #log("Assuming last 2 peaks are relevant: {1}x and {2}x...".format(peaks[1], peaks[2]))
            valid_peaks = peaks[-2:]
            limits = auto_find_limits(valid_peaks, thresholds)
            AUC_low = sum(freqs[: limits["low"]])
            AUC_diplo = sum(freqs[limits["low"] : limits["dip"]])
            AUC_haplo = sum(freqs[limits["dip"] : limits["high"]])
            AUC_high = sum(freqs[limits["high"]:])
            haplotigs_peak_ratio = round(freqs[min(valid_peaks)]/freqs[max(valid_peaks)], 3)


    log("Scoring assembly...")
    AUC_ratio = 1 - (AUC_diplo / AUC_haplo)
    Haploidy = AUC_haplo / ( (AUC_diplo)/2 + AUC_haplo )
    print("AUC(Haploid): H = {}".format(AUC_haplo))
    print("AUC(Diploid): D = {}".format(AUC_diplo))
    print("Ratio: 1-(D/H) = {}".format(round(AUC_ratio, 3)))
    print("Haploidy: H/(H + (D/2)) = {}".format(round(Haploidy, 3)))

    TSS = 1 - abs(SIZE - (AUC_haplo + AUC_diplo / 2)) / SIZE

    write_stats(outfile, AUC_haplo, AUC_diplo, AUC_ratio, Haploidy)

    if plot:
        log("Outputting plots...")
        plot_metrics(
            outfile,
            smoothed,
            peaks,
            heights,
            limits,
            haplotigs_peak_ratio,
            AUC_ratio,
            TSS,
            AUC_low,
            AUC_diplo,
            AUC_haplo,
            AUC_high,
            skip_smooth,
        )

    log("Finished!")
    sys.exit(0)

def write_stats(outname: str, AUC_haplo: float, AUC_diplo: float, AUC_ratio: float,
                haploidy: float
               ):
    log("Outputting stats...")

    f = open(outname + ".stats.txt", "w")
    f.write("AUC(Haploid) = {}\n".format(AUC_haplo))
    f.write("AUC(Diploid) = {}\n".format(AUC_diplo))
    f.write("Ratio = {}\n".format(AUC_ratio))
    f.write("Haploidy = {}\n".format(haploidy))
    f.close()
