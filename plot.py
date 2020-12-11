#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import numpy as np
import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt


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
        "\nHaplotig ratio = {}\nHaploidy = {}\nTSS = {}".format(
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
