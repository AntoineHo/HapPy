#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Requirements:
- sambamba
- scipy
- pandas
- numpy
-
"""

# General modules
import os
import sys
import datetime
from time import localtime, strftime
import argparse

# Subprocessing modules
import subprocess

# from multiprocessing import Pool, TimeoutError

# Signal analysis
import math
import pandas as pd
import numpy as np
from scipy.signal import savgol_filter  # Smoothing
from scipy.signal import find_peaks  # Finding peaks
from scipy.signal import peak_widths

# from scipy.stats import norm
# from scipy.stats import spearmanr
# from scipy.stats import gaussian_kde
# from scipy.optimize import curve_fit

# from Bio import SeqIO # Need BIOPYTHON SEQ/IO
# from pysam import VariantFile # Need Pysam

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt

# from matplotlib.patches import Patch
# from matplotlib.lines import Line2D

from depth import *

#    ______         _     _                       _
#   |  ____|       | |   (_)                     | |
#   | |__     ___  | |_   _   _ __ ___     __ _  | |_    ___
#   |  __|   / __| | __| | | | '_ ` _ \   / _` | | __|  / _ \
#   | |____  \__ \ | |_  | | | | | | | | | (_| | | |_  |  __/
#   |______| |___/  \__| |_| |_| |_| |_|  \__,_|  \__|  \___|
#
#


def estimate_haploidy(args):
    """Finds peaks and modality, then computes scores of haploidy"""
    # Get histogram file
    HIST = check_files(args.HIST)[0]
    # Get size estimation
    is_valid, SIZE = size_from_string(args.SIZE[0])
    if not is_valid:
        log("Error in the size estimation provided...")
        sys.exit(1)

    # Output directory
    outdir = os.path.abspath(os.path.join(os.getcwd(), args.OUT[0]))

    # Get other arguments
    dc_args = {
        "max_contam": args.max_contaminant[0],
        "max_diploid": args.max_diploid[0],
        "min_peak": args.min_peak[0],
        "plot": args.plot,
    }

    _debug = args.debug

    print("# Hap.py estimate")
    print("Coverage histogram:\t{}\n".format(HIST))
    print("Output directory:\t{}\n".format(outdir))
    print("Other arguments: " + str(dc_args))
    print(
        "===============================================================================\n"
    )

    if os.path.isdir(outdir):
        print("WARNING: Output directory already exists! {}".format(outdir))
    else:
        os.makedirs(outdir)  # Create directory following path

    # Read coverage histogram
    log("Reading histogram!")
    freqs = [int(line.split()[1]) for line in open(HIST, "r")][1:-1]

    # Apply Smoothing
    log("Analysing curve!")
    smoothed = [s for s in savgol_filter(freqs, 41, 3) if s >= 0]
    peaks, props = find_peaks(smoothed, height=dc_args["min_peak"])
    heights = props["peak_heights"]
    widths = peak_widths(smoothed, peaks)[0]  # Get peak widths

    if _debug:
        print("freqs:", freqs)
        print("Smoothed:", smoothed)
        print("Peaks:", peaks)
        print("Heights:", heights)
        print("Widths:", widths)

        fig, ax = plt.subplots(figsize=(10, 10))
        # Main distribution
        ax.plot(smoothed, color="k", lw=1.3, zorder=5)
        ax.fill_between(
            np.arange(len(smoothed)),
            smoothed,
            np.zeros(len(smoothed)),
            color="red",
            lw=0,
            alpha=0.35,
            zorder=5,
            label="Smoothed distribution",
        )
        # Vertical lines
        ax.vlines(
            peaks,
            min(smoothed),
            1.06 * max(smoothed),
            label="Peaks found",
            linewidth=1.0,
            zorder=7,
        )
        # Peaks found
        ax.scatter(
            peaks,
            heights,
            marker="x",
            s=100,
            color="b",
            label="Local maximum",
            zorder=9,
        )
        # Formatting plot
        ax.set_xlim(-4, len(smoothed) + 4)
        ax.set_ylim(0, 1.05 * max(smoothed))
        ax.set_xlabel("Coverage", fontsize=15)
        ax.set_ylabel("Frequency", fontsize=15)
        ax.legend()
        ax.locator_params(nbins=40)
        #
        fig.savefig(os.path.join(outdir, "debugplot.png"))
        plt.close(fig)

    params, cov = None, None
    peak_ratios = None
    limits = {}

    if len(peaks) > 3:  # In case 3+ peaks:
        peaks, heights, widths = check_peaks(
            peaks, props, widths, len(smoothed), dc_args
        )

    if len(peaks) == 3:  # In case 3 peaks : 1 in each category
        log("Found 3 peaks at: {}x, {}x and {}x".format(peaks[0], peaks[1], peaks[2]))
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
        log("Found 2 peaks at: {}x and {}x".format(peaks[0], peaks[1]))
        # Peaks could be :
        # CONTA + HAPLO
        # CONTA + DIPLO
        # HAPLO + DIPLO
        haplotigs_peak_ratio = round(heights[-2] / heights[-1], 3)

        if (
            peaks[0] < dc_args["max_contam"]
        ):  # In case found contaminant <-- 2 peaks == contams and higher or lower
            limits["Contaminants"] = math.ceil(peaks[0] + widths[0])

            if peaks[1] < dc_args["max_diploid"]:  # CONTA + DIPLO
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

        if peaks[0] < dc_args["max_contam"]:  # CONTAMINANT ASM
            limits["Contaminants"] = math.ceil(peaks[0] + widths[0])
            limits["Diploid"] = math.ceil(peaks[0] + widths[0])
        elif peaks[0] < dc_args["max_diploid"]:  # DIPLOID ASM
            limits["Contaminants"] = math.ceil(peaks[0] - widths[0])
            limits["Diploid"] = math.ceil(peaks[0] + widths[0])
        else:  # HAPLOID ASM
            limits["Contaminants"] = 0
            limits["Diploid"] = math.ceil(peaks[0] - widths[0])

    else:
        log("No peaks found...")
        log("Finished!")
        sys.exit(0)

    AUC_conta = sum(smoothed[: limits["Contaminants"]])
    AUC_diplo = sum(smoothed[limits["Contaminants"] : limits["Diploid"]])
    AUC_haplo = sum(smoothed[limits["Diploid"] :])

    log("Scoring assembly...")
    AUC_ratio = 1 - (AUC_diplo / AUC_haplo)
    print("AUC(Haploid) (= H) = {}".format(AUC_haplo))
    print("AUC(Diploid) (= D) = {}".format(AUC_diplo))
    print("AUC ratio (1 - D/H) = {}".format(round(AUC_ratio, 3)))

    print("AUC(Haploid) + AUC(Diploid) / 2 = {}".format(AUC_haplo + AUC_diplo / 2))
    TSS = 1 - abs(SIZE - (AUC_haplo + AUC_diplo / 2)) / SIZE
    print("Total Size Score = {}".format(round(TSS, 3)))

    log("Outputting text files...")
    input_basename = os.path.basename(HIST)
    out_file = os.path.join(
        outdir, input_basename + ".stats.txt"
    )  # Create output text filepath

    f = open(out_file, "w")
    f.write("AUC(Haploid) = {}\n".format(AUC_haplo))
    f.write("AUC(Diploid) = {}\n".format(AUC_diplo))
    f.write("AUC ratio (1 - D/H) = {}\n".format(AUC_ratio))
    f.write("AUC(Haploid) + AUC(Diploid)/2 = {}\n".format(AUC_haplo + AUC_diplo / 2))
    f.write("Total Size Score = {}\n".format(TSS))
    f.close()

    if dc_args["plot"]:
        log("Outputting plots...")
        plotfile = os.path.join(
            outdir, input_basename + ".plot.png"
        )  # Create plot filename

        fig, ax = plt.subplots(
            ncols=2,
            nrows=1,
            figsize=(15, 10),
            gridspec_kw={"width_ratios": [0.8, 0.2], "wspace": 0.1},
        )
        # Main distribution
        ax[0].plot(smoothed, color="k", lw=1.3, zorder=5)
        ax[0].fill_between(
            np.arange(len(smoothed)),
            smoothed,
            np.zeros(len(smoothed)),
            color="red",
            lw=0,
            alpha=0.35,
            zorder=5,
            label="Smoothed distribution",
        )
        # Vertical lines
        ax[0].vlines(
            peaks, 0, 1.06 * max(smoothed), label="Peaks found", linewidth=1.0, zorder=7
        )
        for k, v in limits.items():
            col = "r" if k == "Contaminants" else "g"
            ax[0].vlines(
                [v],
                0,
                1.06 * max(smoothed),
                label="Limit for " + k,
                linewidth=1.0,
                color=col,
                zorder=7,
            )
        # Peaks found
        ax[0].scatter(
            peaks,
            heights,
            marker="x",
            s=100,
            color="b",
            label="Local maximum",
            zorder=9,
        )

        # Formatting plot
        ax[0].set_xlim(-4, len(smoothed) + 4)
        ax[0].set_ylim(0, 1.05 * max(smoothed))
        ax[0].set_title(
            "\nHaplotig ratio = {}\nAUC ratio = {}\nTSS = {}".format(
                haplotigs_peak_ratio, round(AUC_ratio, 3), round(TSS, 3)
            ),
            fontsize=16,
        )
        ax[0].set_xlabel("Coverage", fontsize=15)
        ax[0].set_ylabel("Frequency", fontsize=15)
        ax[0].legend()
        ax[0].locator_params(nbins=40)

        # AUC
        ax[1].bar(
            np.arange(3),
            [AUC_conta, AUC_diplo, AUC_haplo],
            color=["darkgray", "violet", "red"],
            edgecolor="k",
            zorder=10,
        )
        ax[1].set_xticks(np.arange(3))
        ax[1].set_xticklabels(
            ["Contam", "Diploid", "Haploid"], rotation="vertical", fontsize=13
        )
        ax[1].set_ylabel("Area Under Curve", fontsize=15)
        ax[1].yaxis.tick_right()
        for axi in ax:
            axi.yaxis.set_tick_params(labelsize=12)
            axi.xaxis.set_tick_params(labelsize=12)
            axi.set_facecolor("whitesmoke")
            axi.yaxis.grid(True, zorder=-1, linewidth=0.5, alpha=0.5)

        fig.savefig(plotfile)
        plt.close(fig)

    log("Finished!")
    sys.exit(0)


def check_peaks(peaks, props, widths, maximum_cov, dc_args):
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
    for pos, height in zip(peaks, props["peak_heights"]):
        if pos < dc_args["max_contam"]:
            if contaminant_peak == None:
                contaminant_peak = pos
                contaminant_height = height
            elif height > contaminant_height:
                contaminant_peak = pos
                contaminant_height = height
            else:
                continue
        elif pos >= dc_args["max_contam"] and pos < dc_args["max_diploid"]:
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
        abs(dc_args["max_contam"] - contaminant_peak)
    )  # floor of absolute distance between contaminant boundary and highest contaminant peak
    diploid_width = math.floor(
        abs(dc_args["max_diploid"] - diploid_peak)
    )  # floor of absolute distance between diploid boundary and highest diploid peak
    haploid_width = math.floor(
        abs(maximum_cov - diploid_peak)
    )  # floor of absolute distance between haploid peak and max coverage ( == len(x) == number of bins ) from the histogram

    """
    # In case number of peaks lower than max contam is > 1
    if len([peak for peak in peaks if peak < dc_args["max_contam"]]) > 1 :
        # INCREASE THRESHOLD to detect peaks
        log("WARNING: More than 1 contaminant peak is found. Trying to use only largest found")

        print("HELP: Try running with --debug flag and check the histogram curve.")
        print("HELP: If you estimate that there should be no contaminant peaks detected then increase the min_peak (-mp) threshold (currently {}).".format(dc_args["min_peak"]))
        print("HELP: If you estimate that there should be a contaminant peak but it is not considered contaminant, then increase the max_contaminant (-mc) optional argument value (currently {}).".format(dc_args["max_contam"]))
        log("Exiting...")
        sys.exit(1)
    elif len([peak for peak in peaks if (peak >= dc_args["max_contam"] and peak < dc_args["max_diploid"])]) > 1 : # If number of diploid peaks is > 1
        log("WARNING: More than 1 diploid peak is found")
        print("HELP: Try running with --debug flag and check the histogram curve.")
        print("HELP: If you estimate that there should be no diploid peaks (or only one) then maybe one is a contaminant or an haploid peak try modifying the max_contaminants and max_diploid arguments (currently {} and {}).".format(dc_args["max_contam"], dc_args["max_diploid"]))
        #print("NOTE: If you estimate that there should be a contaminant peak but it is not considered contaminant, then increase the max_contaminant (-mc) optional argument value (currently {}).".format(dc_args["max_contaminant"]))
        log("Exiting...")
        sys.exit(1)
    else :
        log("WARNING: More than 1 diploid peak is found")
        print("HELP: Try running with --debug flag and check the histogram curve.")
        log("Exiting...")
        sys.exit(1)
    """

    newpeaks = [contaminant_peak, diploid_peak, haploid_peak]
    newheights = [contaminant_height, diploid_height, haploid_height]
    newwidths = [contaminant_width, diploid_width, haploid_width]
    return newpeaks, newheights, newwidths


#    _____    _           _
#   |  __ \  | |         | |
#   | |__) | | |   ___   | |_
#   |  ___/  | |  / _ \  | __|
#   | |      | | | (_) | | |_
#   |_|      |_|  \___/   \__|
#
#


def do_plot():
    pass


#    _____                                  _
#   |  __ \                                | |
#   | |__) |   ___   _ __     ___    _ __  | |_
#   |  _  /   / _ \ | '_ \   / _ \  | '__| | __|
#   | | \ \  |  __/ | |_) | | (_) | | |    | |_
#   |_|  \_\  \___| | .__/   \___/  |_|     \__|
#                   | |
#                   |_|


def build_report():
    pass


#     ____    _______   _    _   ______   _____     _____
#    / __ \  |__   __| | |  | | |  ____| |  __ \   / ____|
#   | |  | |    | |    | |__| | | |__    | |__) | | (___
#   | |  | |    | |    |  __  | |  __|   |  _  /   \___ \
#   | |__| |    | |    | |  | | | |____  | | \ \   ____) |
#    \____/     |_|    |_|  |_| |______| |_|  \_\ |_____/
#
#


def size_from_string(string):
    size_multiplier = {"K": 1000, "M": 1000000, "G": 1000000000}
    if string[-1].upper() in ["K", "M", "G"]:
        if string[:-1].isdigit():
            return True, int(string[:-1]) * size_multiplier[string[-1]]
        else:  # ERROR
            return False, None
    else:  # If no multiplier found
        if string.isdigit():
            return True, int(string)
        else:  # ERROR
            return False, None


def log(string):
    """Prints a log"""
    print("\n{}: {}".format(strftime("%Y-%m-%d %H:%M:%S", localtime()), string))


def check_dirs(dirs):
    """Returns absolute paths and raise exception if dir does not exist"""
    absdirs = []
    for d in dirs:
        if not os.path.isdir(d):
            raise Exception("ERROR: {} is not found!".format(d))
        else:
            absdirs.append(os.path.abspath(d))
    return absdirs


def check_files(files):
    """Returns absolute file paths and raise exception if file does not exist"""
    absfiles = []
    for file in files:
        if not os.path.isfile(file):
            raise Exception("ERROR: {} is not found!".format(file))
        else:
            absfiles.append(os.path.abspath(file))
    return absfiles


def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return False


def list_str(v):
    return v.split(",")


def run(cmd):
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    proc.communicate()
    return proc.returncode


#               _____     _____   _    _   __  __   ______   _   _   _______    _____
#       /\     |  __ \   / ____| | |  | | |  \/  | |  ____| | \ | | |__   __|  / ____|
#      /  \    | |__) | | |  __  | |  | | | \  / | | |__    |  \| |    | |    | (___
#     / /\ \   |  _  /  | | |_ | | |  | | | |\/| | |  __|   | . ` |    | |     \___ \
#    / ____ \  | | \ \  | |__| | | |__| | | |  | | | |____  | |\  |    | |     ____) |
#   /_/    \_\ |_|  \_\  \_____|  \____/  |_|  |_| |______| |_| \_|    |_|    |_____/
#
#


def main():
    """Pipeline check"""
    if not which("sambamba"):
        log("Error: Sambamba is not found!")
        sys.exit(1)

    """Argument parser"""
    parser = argparse.ArgumentParser(
        description="Estimate assembly haploidy based on base depth of coverage histogram."
    )
    subparsers = parser.add_subparsers(required=True, dest="depth")

    # Get coverage of reads
    depth = subparsers.add_parser(
        "depth",
        help="Obtain a depth of coverage histogram for the assembly using sambamba depth base.",
    )
    depth.add_argument(
        "BAM",
        nargs=1,
        type=str,
        help="<STRING> A path to the reads alignment on the assembly (.bam file).",
    )
    depth.add_argument(
        "OUT",
        nargs=1,
        type=str,
        help="<STRING> An output path for the coverage files (.cov and .hist files).",
    )
    depth.add_argument(
        "-t",
        "--threads",
        nargs=1,
        type=int,
        default=[4],
        required=False,
        help="<INT> Number of threads for sambamba. Default: %(default)s.",
    )
    depth.set_defaults(func=get_depth_hist)
    # depth.add_argument('REF',nargs=1,type=str,help="<STRING> A path to the assembly sequence (.fasta file).")

    # Estimate Haploidy and Total Size Score
    estimate = subparsers.add_parser(
        "estimate", help="Computes haploidy score based on the coverage distribution."
    )
    # Positionals
    estimate.add_argument(
        "HIST",
        nargs=1,
        type=str,
        help="<STRING> A path to the histogram output of the `Hap.py depth` command (.hist file).",
    )
    estimate.add_argument(
        "OUT", nargs=1, type=str, help="<STRING> A path for the output directory."
    )
    estimate.add_argument(
        "SIZE",
        nargs=1,
        type=str,
        help="<STRING> An expected assembly size (in bp) to compute the Total Size Score. Valid multipliers are (K, M, G) e.g.: 10K = 10000.",
    )
    # Optionals
    estimate.add_argument(
        "-mc",
        "--max-contaminant",
        nargs=1,
        type=int,
        default=[35],
        required=False,
        help="<INT> Maximum coverage of contaminants. Default: %(default)s",
    )
    estimate.add_argument(
        "-md",
        "--max-diploid",
        nargs=1,
        type=int,
        default=[120],
        required=False,
        help="<INT> Maximum coverage of the diploid peak. Default: %(default)s",
    )
    estimate.add_argument(
        "-mp",
        "--min-peak",
        nargs=1,
        type=int,
        default=[150000],
        required=False,
        help="<INT> Minimum peak height. Default: %(default)s",
    )
    # Flags
    estimate.add_argument(
        "-p",
        "--plot",
        dest="plot",
        action="store_true",
        help="Output plots. Default: %(default)s",
    )
    estimate.add_argument(
        "-d", "--debug", dest="debug", action="store_true", help=argparse.SUPPRESS
    )
    estimate.set_defaults(func=estimate_haploidy, debug=False)
    # estimate.add_argument('FASTA',nargs=1,type=str,help="<STRING> A path to the assembly sequence (.fasta file).")
    # estimate.add_argument('-s', '--sample',nargs=1,type=str,default=['unknown'],help="<STRING> Sample name (only for output file names). Default: %(default)s")
    # estimate.add_argument('-o', '--output',nargs=1,type=str,default=['compare_plots'],help="<STRING> Directory name to output. Default: path/to/cur_dir/%(default)s")
    # estimate.add_argument('-ci', '--color-in',nargs=1,type=str,default=['dodgerblue'],help="<STRING> A matplotlib valid color for the inside distribution. Default: %(default)s")
    # estimate.add_argument('-co', '--color-out',nargs=1,type=str,default=['orange'],help="<STRING> A matplotlib valid color for the outside distribution. Default: %(default)s")
    # estimate.add_argument('-hg','--height',nargs=1,type=int,default=[7], required=False, help="<INT> Height of plot to output. Default: %(default)s")
    # estimate.add_argument('-cs','--contig-size',nargs=1,type=int,default=[500000], required=False, help="<INT> Minimum size of a contig to plot. Default: %(default)sbp")
    # estimate.add_argument('-rs','--region-size',nargs=1,type=int,default=[200], required=False, help="<INT> Minimum size of a region to be considered in output plots. Default: %(default)sbp")
    # estimate.add_argument('-MS','--max-size',nargs=1,type=int,default=[500000], required=False, help="<INT> Maximum size of a region to plot the histogram (regions in range [ms, MS]). Default: %(default)s")
    # estimate.add_argument('-ms','--min-size',nargs=1,type=int,default=[200], required=False, help="<INT> Minimum size of a region to plot the histogram (regions in range [ms, MS]). Default: %(default)s")
    # estimate.add_argument('-bn','--bin-number',nargs=1,type=int,default=[100], required=False, help="<INT> Number of bins in histogram in range [ms, MS]. Default: %(default)s")
    # estimate.add_argument('-MH','--max-het',nargs=1,type=int,default=[10], required=False, help="<INT> Maximum heterozygosity range of a region to plot the histogram (in range [mh, MH]). Default: %(default)s")
    # estimate.add_argument('-mh','--min-het',nargs=1,type=int,default=[0], required=False, help="<INT> Minimum heterozygosity range of a region to plot the histogram (in range [mh, MH]). Default: %(default)s")
    # estimate.add_argument('-bhn','--bin-het-number',nargs=1,type=int,default=[100], required=False, help="<INT> Number of bins in histogram in range [mh, MH]. Default: %(default)s")

    # Plot shuffled vs observed distributions based on shuffle command results
    # plot = subparsers.add_parser('plot', help="Plot peaks and AUC.")
    # plot.add_argument('CSV',nargs=1,type=str,help="<STRING> A path to a the output of the shuffle command (csv file).")
    # plot.add_argument('FASTA',nargs=1,type=str,help="<STRING> A path to the reference genome (fasta file).")
    # plot.add_argument('-o', '--output',nargs=1,type=str,default=['plots'],help="<STRING> Directory name to output. Default: path/to/cur_dir/%(default)s")
    # plot.add_argument('-co', '--color-observed',nargs=1,type=str,default=['dodgerblue'],help="<STRING> A matplotlib valid color for the observed distribution. Default: %(default)s")
    # plot.add_argument('-cs', '--color-shuffled',nargs=1,type=str,default=['orange'],help="<STRING> A matplotlib valid color for the shuffled distribution. Default: %(default)s")
    # plot.add_argument('-wd','--width',nargs=1,type=int,default=[14], required=False, help="<INT> Width of plot to output. Default: %(default)s")
    # plot.add_argument('-hg','--height',nargs=1,type=int,default=[7], required=False, help="<INT> Height of plot to output. Default: %(default)s")
    # plot.add_argument('-s','--size',nargs=1,type=int,default=[500000], required=False, help="<INT> Minimum size of a contig to plot. Default: %(default)sbp")
    # plot.add_argument('-b','--bins',nargs=1,type=int,default=[501], required=False, help="<INT> Number of bins (windows) to make on each chromosome. Default: %(default)sbp")
    # plot.set_defaults(func=do_plot)

    # Creates a report for multiple assemblies
    # report = subparsers.add_parser('report', help="Compare heterozygosity of regions in and out of a bed file.")
    # report.add_argument('VCF',nargs=1,type=str,help="<STRING> A path to a the short variants calls (vcf file).")
    # report.add_argument('FASTA',nargs=1,type=str,help="<STRING> A path to the reference genome (fasta file).")
    # report.add_argument('BED',nargs=1,type=str,help="<STRING> A path to the regions locations (sorted and merged bed file).")
    # report.add_argument('-s', '--sample',nargs=1,type=str,default=['ancestor'],help="<STRING> Sample name to read in VCF. Default: %(default)s")
    # report.add_argument('-o', '--output',nargs=1,type=str,default=['plots'],help="<STRING> Directory name to output. Default: path/to/cur_dir/%(default)s")
    # report.add_argument('-ci', '--color-in',nargs=1,type=str,default=['dodgerblue'],help="<STRING> A matplotlib valid color for the inside distribution. Default: %(default)s")
    # report.add_argument('-co', '--color-out',nargs=1,type=str,default=['orange'],help="<STRING> A matplotlib valid color for the outside distribution. Default: %(default)s")
    # report.add_argument('-wd','--width',nargs=1,type=int,default=[14], required=False, help="<INT> Width of plot to output. Default: %(default)s")
    # report.add_argument('-hg','--height',nargs=1,type=int,default=[7], required=False, help="<INT> Height of plot to output. Default: %(default)s")
    # report.add_argument('-cs','--contig-size',nargs=1,type=int,default=[500000], required=False, help="<INT> Minimum size of a contig to plot. Default: %(default)sbp")
    # report.add_argument('-rs','--region-size',nargs=1,type=int,default=[200], required=False, help="<INT> Minimum size of a region to be considered in output plots. Default: %(default)sbp")
    # report.add_argument('-MS','--max-size',nargs=1,type=int,default=[500000], required=False, help="<INT> Maximum size of a region to plot the histogram (regions in range [ms, MS]). Default: %(default)s")
    # report.add_argument('-ms','--min-size',nargs=1,type=int,default=[200], required=False, help="<INT> Minimum size of a region to plot the histogram (regions in range [ms, MS]). Default: %(default)s")
    # report.add_argument('-bn','--bin-number',nargs=1,type=int,default=[100], required=False, help="<INT> Number of bins in histogram in range [ms, MS]. Default: %(default)s")
    # report.add_argument('-MH','--max-het',nargs=1,type=int,default=[10], required=False, help="<INT> Maximum heterozygosity range of a region to plot the histogram (in range [mh, MH]). Default: %(default)s")
    # report.add_argument('-mh','--min-het',nargs=1,type=int,default=[0], required=False, help="<INT> Minimum heterozygosity range of a region to plot the histogram (in range [mh, MH]). Default: %(default)s")
    # report.add_argument('-bhn','--bin-het-number',nargs=1,type=int,default=[100], required=False, help="<INT> Number of bins in histogram in range [mh, MH]. Default: %(default)s")
    # report.set_defaults(func=build_report)

    args = parser.parse_args()
    args.func(args)

    sys.exit(0)


if __name__ == "__main__":
    main()

