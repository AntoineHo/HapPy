#!/usr/bin/python
# -*- coding: utf-8 -*-

# General
import os

# Stats & representation
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

# Happy
from happy.utils import *

# DEBUGGING PLOT
def debug_smooth_histogram(freqs, smoothed, peaks, heights, widths, outdir):
    """Plot an histogram to adjust options"""
    # NOTE useless here
    #print("freqs:", freqs)
    #print("Smoothed:", smoothed)
    #print("Peaks:", peaks)
    #print("Heights:", heights)
    #print("Widths:", widths)

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
        peaks, heights, marker="x", s=100, color="b", label="Local maximum", zorder=9
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

# AUTO ESTIMATE MODULE
def plot_model(x, smoothed, curve_function, popt, peak_opts,
               expected_peaks, actual_peaks, peak_areas,
               peak_multipliers, outdir, outname,
               ) :

    fig, axs = plt.subplots(ncols=3, nrows=1, figsize=(15,5), gridspec_kw={"width_ratios":[5,5,5]})
    ax = axs[0]
    ax.plot(x, smoothed, 'b-', label='Observed')
    ax.plot(x, curve_function(x, *popt), 'r-', label='Model')

    for peak, opts in peak_opts.items() :
        curve = gauss(x, *opts)
        ax.plot(x, curve, linewidth=0.5, color="k")
        ax.fill_between(x, curve, alpha=0.5, label="{}X".format(peak))

    ax.vlines(expected_peaks, ymin=0, ymax=max(smoothed), linestyle="-", color="r", linewidth=0.75)
    ax.vlines(actual_peaks, ymin=0, ymax=max(smoothed), linestyle="--", color="k", linewidth=0.75)
    ax.legend()

    ax = axs[1]
    ax.plot(x, smoothed, 'b-', label='Observed')
    ax.plot(x, curve_function(x, *popt), 'r-', label='Model')
    ax.fill_between(x, curve_function(x, *popt), alpha=0.5)
    ax.legend()

    ax = axs[2]
    other_area = sum(smoothed) - sum(peak_areas.values())
    ax.bar(x=[0], height=other_area, edgecolor="k")
    xticks = [0]
    xticklabels = ["Others"]
    i = 1
    for peak, area in peak_areas.items() :
        ax.bar([i, i+1], height=[area, area*peak_multipliers[peak]], edgecolor="k")
        xticks += [i, i+1]
        xticklabels += [str(peak)+"X", "Divided"]
        i += 2

    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels, rotation="vertical")

    for ax in axs :
        ax.yaxis.set_tick_params(labelsize=12)
        ax.xaxis.set_tick_params(labelsize=12)
        ax.set_facecolor('whitesmoke')
        ax.yaxis.grid(True, zorder=-1, linewidth=0.5, alpha=0.5)

    plt.tight_layout()
    fig.savefig(os.path.join(outdir, outname+".png"))
    plt.close(fig)















# ESTIMATE MODULE
def plot_curves(
    combination,
    actual_curve,
    expected_curve,
    ploidy,
    peak_limits,
    determined_peaks_ordered,
    outname,
) :

    peak_limits_plot = [v for k, v in peak_limits.items()]
    peak_limits_plot = list(set(sum(peak_limits_plot, ())))

    fig, ax = plt.subplots(figsize=(10,5))
    x = np.arange(len(actual_curve))
    ax.plot(x, actual_curve, color="orange", alpha=0.85, label="Actual curve")
    ax.plot(x, expected_curve, linestyle="--", color="k", linewidth=0.75, label="Determined curve")
    ax.vlines(peak_limits_plot, min(actual_curve), max(actual_curve))
    ax.vlines(determined_peaks_ordered, min(actual_curve), max(actual_curve), color="k", linewidth=0.75, linestyle="--")

    for peak, limits in peak_limits.items() :
        x = np.arange(len(actual_curve))
        regions = [(i >= limits[0]) & (i < limits[1]) for i in x]

        if peak in determined_peaks_ordered :
            peakploidy = ploidy - determined_peaks_ordered.index(peak)
            label = "Peak {}X".format(peakploidy)
        else :
            label = "Contaminant"
        ax.fill_between(x, actual_curve, where=regions, alpha=0.5, label=label)

    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.15), fancybox=True, shadow=False, ncol=5)

    title_string = "Organism peaks: ({})\n".format(",".join(str(i) for i in combination["Peaks"]))
    title_string += "Contaminants: ({})\n".format(",".join(str(i) for i in combination["Contam"]))
    title_string += "Determined peaks: ({})".format(",".join(str(i) for i in determined_peaks_ordered))
    ax.set_title(title_string)

    plt.show()


def plot_metrics(
    outname,
    smoothed_freq,
    peaks,
    heights,
    limits,
    haplotigs_peak_ratio,
    AUC_ratio,
    TSS,
    AUC_conta,
    AUC_diplo,
    AUC_haplo,
):

    plotfile = os.path.join(outname + ".plot.png")  # Create plot filename

    fig, ax = plt.subplots(
        ncols=2,
        nrows=1,
        figsize=(15, 10),
        gridspec_kw={"width_ratios": [0.8, 0.2], "wspace": 0.1},
    )
    # Main distribution
    ax[0].plot(smoothed_freq, color="k", lw=1.3, zorder=5)
    ax[0].fill_between(
        np.arange(len(smoothed_freq)),
        smoothed_freq,
        np.zeros(len(smoothed_freq)),
        color="red",
        lw=0,
        alpha=0.35,
        zorder=5,
        label="Smoothed distribution",
    )
    # Vertical lines
    ax[0].vlines(
        peaks,
        0,
        1.06 * max(smoothed_freq),
        label="Peaks found",
        linewidth=1.0,
        zorder=7,
    )
    for k, v in limits.items():
        col = "r" if k == "Contaminants" else "g"
        ax[0].vlines(
            [v],
            0,
            1.06 * max(smoothed_freq),
            label="Limit for " + k,
            linewidth=1.0,
            color=col,
            zorder=7,
        )
    # Peaks found
    ax[0].scatter(
        peaks, heights, marker="x", s=100, color="b", label="Local maximum", zorder=9
    )

    # Formatting plot
    ax[0].set_xlim(-4, len(smoothed_freq) + 4)
    ax[0].set_ylim(0, 1.05 * max(smoothed_freq))
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
