#!/usr/bin/python
# -*- coding: utf-8 -*-

from time import localtime, strftime
import os, subprocess
import numpy as np


def log(string: str):
    """Prints a log"""
    print("\n{}: {}".format(strftime("%Y-%m-%d %H:%M:%S", localtime()), string))


def run(cmd):
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    proc.communicate()
    return proc.returncode


def check_files(file):
    """Returns absolute file paths and raise exception if file does not exist"""
    if not os.path.isfile(file):
        raise Exception("ERROR: {0} is not found!".format(file))
    return file


def size_from_string(string: str) -> int:
    size_multiplier = {"K": 1000, "M": 1000000, "G": 1000000000}
    if string[-1].upper() in ["K", "M", "G"] and string[:-1].isdigit():
        return int(string[:-1]) * size_multiplier[string[-1]]
    elif string.isdigit():
        return int(string)
    else:  # ERROR
        raise Exception("Size argument is not a valid number.")


def debug_smooth_histogram(freqs, smoothed, peaks, heights, widths):
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
