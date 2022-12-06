#!/usr/bin/python
# -*- coding: utf-8 -*-

import argparse

def execute_coverage(args) :
    """Coverage histogram command"""
    try :
        import happy.coverage as happyc
    except :
        import coverage as happyc

    happyc.get_cov_hist(
        args.MAP[0], args.threads[0], args.outdir[0], args.diploid, args.samtools[0],
    )

def execute_estimate(args) :
    """Estimate command"""
    try :
        import happy.estimate as happye
    except :
        import estimate as happye

    happye.estimate_haploidy(
        args.COV[0], args.size[0], args.outstats[0],
        args.limit_low[0], args.limit_diploid[0], args.limit_high[0],
        args.min_peak[0], args.prominence[0],
        args.valley_height[0], args.valley_prominence[0],
        args.window[0],
        plot=args.plot,
        debug=args.debug,
        skip_smooth=args.no_smooth,
    )

def execute_autoest(args) :
    """Auto estimate command"""
    try :
        import happy.autoestimate as happyae
    except :
        import autoestimate as happyae

    happyae.auto_estimate_haploidy(
        args.COV[0],
        args.outstats[0],
        args.size[0],
        args.min_peak[0],
        args.prominence[0],
        args.window[0],
        args.score[0],
        args.plot,
        args.debug,
    )


def main():

    """Argument parser"""
    parser = argparse.ArgumentParser(description='Estimate assembly haploidy based on base depth of coverage histogram.')
    subparsers = parser.add_subparsers(required=True, dest="coverage || estimate || autoest")

    cov = subparsers.add_parser('coverage', help="Compute coverage histogram for mapping file.")
    cov.add_argument("MAP",             nargs=1, type=str,                       help="<FILE> Sorted BAM file after mapping reads to the assembly.")
    cov.add_argument("-t", "--threads", nargs=1, type=int, default=[1],          help="<INT> Number of parallel threads allocated for sambamba. Default: %(default)s")
    cov.add_argument("-d", "--outdir",  nargs=1, type=str, default=["out"],      help="<DIR> Path where the .cov and .hist files are written. Default: %(default)s")
    cov.add_argument("--diploid",       action="store_true",                     help="Allows for multimapping reads in sambamba filters. Default: %(default)s")
    cov.add_argument("--samtools",      nargs=1, type=str, default=["samtools"], help="<PATH> When using the --diploid flag only, which samtools executable to use. Due to a problem in sambamba depth base with -F 'mapping quality >= 0', samtools is required to obtain coverage in diploid assemblies.")
    cov.set_defaults(func=execute_coverage)

    est = subparsers.add_parser('estimate', help="Compute haploidy from coverage histogram.")
    est.add_argument("COV", nargs=1, help="<FILE> Coverage histogram file.")
    est.add_argument("-s", "--size",             nargs=1, type=str, required=True,   help="<STRING> Estimated haploid genome size. (Recognized modifiers: K,M,G).")
    est.add_argument("-o", "--outstats",         nargs=1, type=str, required=True,   help="<FILE> Path where haploidy value is written.")
    est.add_argument("-m", "--min-peak",         nargs=1, default=[15000], type=int, help="<INT> Minimum peak height (see SciPy doc). Default: %(default)s")
    est.add_argument("-p", "--prominence",       nargs=1, default=[10000], type=int, help="<INT> Minimum peak prominence (see SciPy docs). Default: %(default)s")
    est.add_argument('-vh','--valley-height',    nargs=1, default=[15000], type=int, help="<INT> Minimum valley height (abs). Default: %(default)s")
    est.add_argument('-vp','--valley-prominence',nargs=1, default=[10000], type=int, help="<INT> Minimum valley prominence (abs). Default: %(default)s")
    est.add_argument('-w','--window',            nargs=1, default=[41],    type=int, help="<INT> Window length to use in savgol filter. Default: %(default)s")
    est.add_argument("-ll", "--limit-low",       nargs=1, default=[None],  type=int, help="<INT> Lower threshold of coverage (lower than diploid peak). Default: auto.")
    est.add_argument("-ld",  "--limit-diploid",  nargs=1, default=[None],  type=int, help="<INT> Middle threshold of coverage (between diploid and haploid peaks). Default: auto.")
    est.add_argument("-lh", "--limit-high",      nargs=1, default=[None],  type=int, help="<INT> Upper threshold of coverage (higher than haploid peak). Default: auto.")
    est.add_argument("--plot",                   action="store_true",                help="Generate histogram plot. Default: %(default)s")
    est.add_argument("--debug",                  action="store_true",                help="Generate more informative plots to debug. Default: %(default)s")
    est.add_argument("--no-smooth",              action="store_true",                help="Skip the smoothing step. Default: %(default)s")
    est.set_defaults(func=execute_estimate)

    auto = subparsers.add_parser('autoest',     help="Detect peaks and computes haploidy metrics from the coverage histogram.")
    auto.add_argument('COV', nargs=1, help="<FILE> Coverage histogram file.")
    auto.add_argument("-s", "--size",       nargs=1,                            help="<INT> Estimated haploid genome size. (Recognized modifiers: K,M,G).")
    auto.add_argument("-o", "--outstats",   nargs=1,                            help="<FILE> Path to file where the metrics values will be written.")
    auto.add_argument("-m", "--min-peak",   nargs=1, default=[15000],type=int,  help="<INT> Minimum peak height. Default: %(default)s")
    auto.add_argument("-p", "--prominence", nargs=1, default=[10000],type=int,  help="<INT> Minimum peak prominence (see SciPy docs). Default: %(default)s")
    auto.add_argument("-w", "--window",     nargs=1, default=[1.5],  type=float,help="<FLOAT> Window size for peak matching modifier. Default: %(default)s")
    auto.add_argument("-sc", "--score",     nargs=1, default=[0.75], type=float,help="<FLOAT> Score threshold for outputting to file. Default: %(default)s")
    auto.add_argument("--plot",             action="store_true",                help="Generate plots. Default: %(default)s")
    auto.add_argument("--debug",            action="store_true",                help="Generate debug histogram plot. Default: %(default)s")
    auto.set_defaults(func=execute_autoest)

    args = parser.parse_args()
    args.func(args)
    sys.exit(0)

if __name__ == "__main__":
    main()
