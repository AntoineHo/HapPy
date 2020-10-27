#!/usr/bin/python
# -*- coding: utf-8 -*-

import os, math
from utils import *
from scipy.signal import savgol_filter  # Smoothing
from scipy.signal import find_peaks  # Finding peaks
from scipy.signal import peak_widths


def estimate_haploidy(
    infile,
    max_cont: int,
    max_dip: int,
    min_peak: int,
    size: int,
    outfile,
    plot=True,
    debug=False,
):
    """Finds peaks and modality, then computes scores of haploidy"""

    # Get histogram file and size estimation
    HIST = check_files(infile)
    SIZE = size_from_string(size)

    # Get other arguments
    dc_args = {
        "max_contam": max_cont,
        "max_diploid": max_dip,
        "min_peak": min_peak,
        "plot": plot,
    }

    print("# Hap.py estimate")
    print(
        "Coverage histogram:\t{0}\nOutput file:\t{1}\nOther arguments:\t{2}\n".format(
            HIST, outfile, dc_args
        )
    )
    print(
        "===============================================================================\n"
    )

    # Read coverage histogram
    log("Reading histogram!")
    freqs = [int(line.split()[1]) for line in open(HIST, "r")][1:-1]

    # Apply Savitzky-Golay filter to smooth coverage histogram
    log("Analysing curve!")
    smoothed = [s for s in savgol_filter(freqs, 41, 3) if s >= 0]
    peaks, props = find_peaks(smoothed)
    heights = [smoothed[i] for i in peaks]
    widths = peak_widths(smoothed, peaks)[0]  # Get peak widths

    if debug:
        debug_smooth_histogram(freqs, smoothed, peaks, heights, widths)

    params, cov = None, None
    peak_ratios = None
    limits = {}

    # if len(peaks) > 3:  # In case 3+ peaks:
    #    peaks, heights, widths = check_peaks(
    #        peaks, props, widths, len(smoothed), dc_args
    #    )

    if len(peaks) == 3:  # In case 3 peaks : 1 in each category
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
