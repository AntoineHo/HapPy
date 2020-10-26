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

