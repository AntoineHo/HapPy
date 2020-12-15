#!/usr/bin/python
# -*- coding: utf-8 -*-
# Based on Rémy Greinhofer (rgreinho) tutorial on subcommands in docopt
# https://github.com/rgreinho/docopt-subcommands-example
# Also based on Cyril Mathey-Doret (cmdoret) commands.py in hicstuff
# https://github.com/koszullab/hicstuff

# General
import sys, os, shutil
import tempfile
from os.path import join, dirname

# Commands
from docopt import docopt

# Happy
import happy.coverage as happyc
import happy.estimate as happye
import happy.autoestimate as happyae


class AbstractCommand:
    """Base class for the commands"""

    def __init__(self, command_args, global_args):
        """Initialize the commands."""
        self.args = docopt(self.__doc__, argv=command_args)
        self.global_args = global_args

    def execute(self):
        """Execute the commands"""
        raise NotImplementedError

    def check_output_path(self, path, force=False):
        """Throws error if the output file exists. Create required file tree otherwise."""
        # Get complete output filename and prevent overwriting unless force is enabled
        if not force and os.path.exists(path):
            raise IOError("Output file already exists. Use --force to overwrite")
        if dirname(path):
            os.makedirs(dirname(path), exist_ok=True)


class Coverage(AbstractCommand):
    """Coverage histogram command
    Compute coverage histogram for mapping file.

    usage:
        coverage [--threads=1] --outdir=DIR <mapping.bam>

    arguments:
        mapping.bam              Sorted BAM file after mapping reads to the assembly.

    options:
        -t, --threads=INT        Number of parallel threads allocated for
                                 sambamba [default: 1].
        -d, --outdir=DIR         Path where the .cov and .hist files are written.
    """

    def execute(self):
        print("Running coverage module.")
        happyc.get_cov_hist(
            self.args["<mapping.bam>"], self.args["--threads"], self.args["--outdir"]
        )


class Estimate(AbstractCommand):
    """Estimate command
    Compute haploidy from coverage histogram.

    usage:
        estimate [--max-contaminant=INT] [--max-diploid=INT] --size=INT --outstats=FILE [--plot] <coverage.hist>

    arguments:
        coverage.hist               Coverage histogram.

    options:
        -C, --max-contaminant=INT   Maximum coverage of contaminants.
        -D, --max-diploid=INT       Maximum coverage of the diploid peak.
        -S, --size=INT              Estimated haploid genome size.
        -O, --outstats=FILE         Path where haploidy value is written.
        -P, --plot                  Generate histogram plot.
    """

    def execute(self):

        happye.estimate_haploidy(
            self.args["<coverage.hist>"],
            self.args["--max-contaminant"],
            self.args["--max-diploid"],
            self.args["--size"],
            self.args["--outstats"],
        )

class Autoest(AbstractCommand):
    """Auto estimate command
    Detect peaks and computes haploidy metrics from the coverage histogram.

    usage:
        autoest [--min-peak=INT] [--prominence=INT] [--window=FLOAT] [--score=FLOAT] [--plot] [--debug] --size=INT --outstats=FILE <coverage.hist>

    arguments:
        coverage.hist               Coverage histogram.

    options:
        -S, --size=STRING           Estimated haploid genome size
                                    (Recognized modifiers: K,M,G).
        -O, --outstats=FILE         Path to file where the metrics values will be written.
        -M, --min-peak=INT          Minimum peak height
                                    [default: 15000].
        -P, --prominence=INT        Minimum peak prominence (see SciPy docs)
                                    [default: 10000].
        -W, --window=FLOAT          Window size for peak matching modifier
                                    [default: 1.5].
        -sc, --score=FLOAT          Score threshold for outputting to file
                                    [default: 0.75].
        -p, --plot                  Generate plots.
        -d, --debug                 Generate debug histogram plot.
    """

    def execute(self):

        happyae.auto_estimate_haploidy(
            self.args["<coverage.hist>"],
            self.args["--outstats"],
            self.args["--size"],
            int(self.args["--min-peak"]),
            int(self.args["--prominence"]),
            float(self.args["--window"]),
            float(self.args["--score"]),
            self.args["--plot"],
            self.args["--debug"],
        )
