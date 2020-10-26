#!/usr/bin/python
# -*- coding: utf-8 -*-
# Based on Rémy Greinhofer (rgreinho) tutorial on subcommands in docopt
# https://github.com/rgreinho/docopt-subcommands-example
# Also based on Cyril Mathey-Doret (cmdoret) commands.py in hicstuff
# https://github.com/koszullab/hicstuff

from docopt import docopt
import sys, os, shutil
import tempfile
from os.path import join, dirname
import coverage as happyc
import estimate as happye


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
        -d, --outdir=DIR           Path where the .cov and .hist files are written.
    """

    def execute(self):
        print("Running coverage module.")
        happyc.get_cov_hist(
            self.args["<mapping.bam>"], self.args["--threads"], self.args["--outdir"]
        )


class Estimate(AbstractCommand):
    """Estimate command
    Compute AUC ratio and TSS from coverage histogram.

    usage:
        estimate --max-contaminant=INT --max-diploid=INT --min-peak=INT --outstats=FILE [--plot] <coverage.hist>
        
    arguments:
        coverage.hist                  Coverage histogram.
        
    options:
        -C, --max-contaminant=INT   Maximum coverage of contaminants.
        -D, --max-diploid=INT       Maximum coverage of the diploid peak.
        -M, --min-peak=INT          Minimum peak height.
        -O, --outstats=FILE            Path where the AUC ratio and TSS values are written.
        -p, --plot                  Generate histogram plot.
    """

    def execute(self):

        happye.estimate_haploidy(
            self.args["<coverage.hist>"],
            self.args["--max-contaminant"],
            self.args["--max-diploid"],
            self.args["--min-peak"],
            self.args["--outstats"],
        )
