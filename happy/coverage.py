#!/usr/bin/python
# -*- coding: utf-8 -*-

# General
import os, sys

# Happy
from happy.utils import *

# Stats
import pandas as pd
import numpy as np


def get_cov_hist(infile, threads: int, outdir):
    """Finds peaks and modality, then computes scores of haploidy"""

    print("# Hap.py coverage")
    print(
        "Input alignment file:\t{0}\nOutput directory:\t{1}\nNumber of threads:\t{2}".format(
            infile, outdir, threads
        )
    )
    print(
        "===============================================================================\n"
    )

    if os.path.isdir(outdir):
        print("WARNING: Output directory already exists! {}".format(outdir))
    else:
        os.makedirs(outdir)  # Create directory following path

    if not which("sambamba") :
        log("ERROR: Sambamba is not found!")
        sys.exit(1)

    # Read coverage histogram
    coverage_output = os.path.join(outdir, os.path.basename(infile) + ".cov")
    if not os.path.isfile(coverage_output):  # In case no coverage file found
        log("Starting sambamba depth...")
        coverage_output = os.path.abspath(coverage_output)
        dc_sambamba = {"BAM": infile, "threads": threads, "out": coverage_output}
        cmd = "sambamba depth base -t {threads} -o {out} --min-coverage=0 --min-base-quality=0 {BAM}"
        cmd = cmd.format(**dc_sambamba)
        sambamba_returncode = run(cmd)
        if sambamba_returncode == 1:
            log(
                "ERROR: sambamba command returned: {0}, a common problem is a missing index (.bai) file...".format(
                    sambamba_returncode
                )
            )
            sys.exit(1)
    else:
        log(
            "SKIP: Existing output coverage file found using it instead of running sambamba!"
        )

    # Read coverage histogram
    output = os.path.join(outdir, os.path.basename(infile) + ".hist")
    if not os.path.isfile(output):  # In case no coverage file found
        log("Reading coverage file...")
        output = os.path.abspath(output)

        # Read coverage
        df = pd.read_csv(coverage_output, sep="\t", usecols=["REF", "POS", "COV"])
        total_bases = len(df)  # ~length of reference assembly
        # Number of bins = number of coverage values between 0 and max found
        hist, bin_edges = np.histogram(df["COV"], bins=range(0, max(df["COV"])))

        total_summarized = 0
        max_coverage_to_keep = 0
        for freq, cov in zip(hist, bin_edges):
            total_summarized += freq
            if total_summarized / total_bases > 0.99:
                max_coverage_to_keep = cov  # Find max coverage to keep
                break

        # Output histogram
        f = open(output, "w")
        for freq, cov in zip(hist, bin_edges):
            if cov == max_coverage_to_keep:
                break
            else:
                f.write("{}\t{}\n".format(cov, freq))
        f.close()

    else:
        log("SKIP: Existing output histogram file found!")

    log("Finished!")
    sys.exit(0)
