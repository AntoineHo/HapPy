#!/usr/bin/python
# -*- coding: utf-8 -*-

# General
import os, sys, math

# Happy
from happy.utils import *
from happy.plot import *

# Stats
from scipy.signal import savgol_filter  # SmoothingAUC
from scipy.signal import find_peaks  # Finding peaks
from scipy.signal import peak_widths

def estimate_haploidy(
    infile, max_cont: int, max_dip: int, size: int, outfile, plot=False, debug=False
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
    log("Reading histogram!")
    freqs = [int(line.split()[1]) for line in open(HIST, "r")][1:-1]

    # Apply Savitzky-Golay filter to smooth coverage histogram
    log("Analysing curve!")
    smoothed = [s for s in savgol_filter(freqs, 41, 3) if s >= 0]
    peaks, props = find_peaks(smoothed, height=15000, prominence=10000)
    heights = [smoothed[i] for i in peaks]
    widths = peak_widths(smoothed, peaks)[0]  # Get peak widths

    if max_cont is None:
        max_cont = peaks[len(peaks) - 1] * 0.30
    else :
        try :
            max_cont = int(max_cont)
        except :
            raise Exception("ERROR: Invalid value of --max-cont!")

    if max_dip is None:
        max_dip = peaks[len(peaks) - 1] * 0.55
    else :
        try :
            max_dip = int(max_dip)
        except :
            raise Exception("ERROR: Invalid value of --max-dip!")

    if debug:
        debug_smooth_histogram(freqs, smoothed, peaks, heights, widths, os.getcwd())

    params, cov = None, None
    peak_ratios = None
    limits = {}

    if len(peaks) > 3:  # In case 3+ peaks:
        peaks, heights, widths = check_peaks(
            peaks, heights, widths, len(smoothed), max_cont, max_dip
        )

    if len(peaks) == 0:
        log("No peak found.")
        sys.exit(1)

    elif len(peaks) == 3:  # In case 3 peaks : 1 in each category
        log(
            "Found 3 peaks at: {0}x, {1}x and {2}x".format(peaks[0], peaks[1], peaks[2])
        )
        haplotigs_peak_ratio = round(heights[-2] / heights[-1], 3)

        # For each peak pos and width
        for n, pos, width in zip(range(3), peaks, widths):
            if n == 0:
                categ = "Contaminants"
            elif n == 1:
                categ = "Diploid"
            elif n == 2:
                break
            # Store limit under the key category of the peak
            limits[categ] = math.ceil(pos + width)

    elif len(peaks) == 2:  # If only 2 peaks
        log("Found 2 peaks at: {0}x and {1}x".format(peaks[0], peaks[1]))
        # Peaks could be :
        # CONTA + HAPLO
        # CONTA + DIPLO
        # HAPLO + DIPLO
        haplotigs_peak_ratio = round(heights[-2] / heights[-1], 3)

        if (
            peaks[0] < max_cont
        ):  # In case found contaminant <-- 2 peaks == contams and higher or lower
            limits["Contaminants"] = math.ceil(peaks[0] + widths[0])

            if peaks[1] < max_dip:  # CONTA + DIPLO
                limits["Diploid"] = math.ceil(
                    peaks[1] - widths[1]
                )  # diploid region is between contams and lower border of haploid peak
            else:  # CONTA + HAPLO
                limits["Diploid"] = math.ceil(
                    peaks[1] + widths[1]
                )  # diploid region is between contams and upper border of diploid peak

        else:  # DIPLO + HAPLO
            limits["Contaminants"] = 0
            limits["Diploid"] = math.ceil(peaks[0] + widths[0])

    elif len(peaks) == 1:
        log("Found 1 peak at: {}x".format(peaks[0]))
        haplotigs_peak_ratio = 0.0

        if peaks[0] < max_cont:  # CONTAMINANT ASM
            limits["Contaminants"] = math.ceil(peaks[0] + widths[0])
            limits["Diploid"] = math.ceil(peaks[0] + widths[0])
        elif peaks[0] < max_dip:  # DIPLOID ASM
            limits["Contaminants"] = math.ceil(peaks[0] - widths[0])
            limits["Diploid"] = math.ceil(peaks[0] + widths[0])
        else:  # HAPLOID ASM
            limits["Contaminants"] = 0
            limits["Diploid"] = math.ceil(peaks[0] - widths[0])

    AUC_conta = sum(smoothed[: limits["Contaminants"]])
    AUC_diplo = sum(smoothed[limits["Contaminants"] : limits["Diploid"]])
    AUC_haplo = sum(smoothed[limits["Diploid"] :])

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
            AUC_conta,
            AUC_diplo,
            AUC_haplo,
        )

    log("Finished!")
    sys.exit(0)


def check_peaks(peaks, heights, widths, maximum_cov: int, max_cont: int, max_dip: int):
    """In case there are more than 3 peaks, find only the 3 highest interest peaks"""
    log(
        "Warning: detected more than 3 peaks at: {}x and {}x".format(
            "x, ".join(str(peak) for peak in peaks[:-1]), peaks[-1]
        )
    )

    contaminant_peak, diploid_peak, haploid_peak = None, None, None
    contaminant_height, diploid_height, haploid_height = None, None, None
    contaminant_width, diploid_width, haploid_width = None, None, None

    # Find highest peaks
    for pos, height in zip(peaks, heights):
        if pos < max_cont :
            if contaminant_peak == None:
                contaminant_peak = pos
                contaminant_height = height
            elif height > contaminant_height:
                contaminant_peak = pos
                contaminant_height = height
            else:
                continue
        elif pos >= max_cont and pos < max_dip:
            if diploid_peak == None:
                diploid_peak = pos
                diploid_height = height
            elif height > diploid_height:
                diploid_peak = pos
                diploid_height = height
            else:
                continue
        else:
            if haploid_peak == None:
                haploid_peak = pos
                haploid_height = height
            elif height > haploid_height:
                haploid_peak = pos
                haploid_height = height
            else:
                continue

    # widths = general boundaries <-- may be a problem sometimes
    contaminant_width = math.floor(
        abs(max_cont - contaminant_peak)
    )  # floor of absolute distance between contaminant boundary and highest contaminant peak
    diploid_width = math.floor(
        abs(max_dip - diploid_peak)
    )  # floor of absolute distance between diploid boundary and highest diploid peak
    haploid_width = math.floor(
        abs(maximum_cov - diploid_peak)
    )  # floor of absolute distance between haploid peak and max coverage ( == len(x) == number of bins ) from the histogram

    newpeaks = [contaminant_peak, diploid_peak, haploid_peak]
    newheights = [contaminant_height, diploid_height, haploid_height]
    newwidths = [contaminant_width, diploid_width, haploid_width]
    return newpeaks, newheights, newwidths

def write_stats(outname: str, AUC_haplo: float, AUC_diplo: float, AUC_ratio: float,
                haploidy: float
               ):
    log("Outputting stats...")

    f = open(outname, "w")
    f.write("AUC(Haploid) = {}\n".format(AUC_haplo))
    f.write("AUC(Diploid) = {}\n".format(AUC_diplo))
    f.write("Ratio = {}\n".format(AUC_ratio))
    f.write("Haploidy = {}".format(haploidy))
    f.close()
