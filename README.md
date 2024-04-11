# HapPy
**Hap**loidy using **Py**thon.

Easy haploidy estimation.

[![DOI](https://zenodo.org/badge/299235590.svg)](https://zenodo.org/badge/latestdoi/299235590)
[![MIT License](https://img.shields.io/github/license/AntoineHo/HapPy)](https://github.com/AntoineHo/HapPy/blob/master/LICENSE)
[![PyPI version](https://badge.fury.io/py/happy-AntoineHo.svg)](https://badge.fury.io/py/happy-AntoineHo)

## 1. General
This tool assesses the haploidy *H* of a given assembly.
*H* is defined as the fraction of the bases of the genome that are in the collapsed peak *C*. This metrics is calculated as *H*=*C*/(*C*+*U*/2), where *C* is the size (area) of the collapsed peak and *U* the size of the uncollapsed peak in the per-base coverage histogram of the assembly.

For more information, see:
  **Overcoming uncollapsed haplotypes in long-read assemblies of non-model organisms**,
  Nadège Guiglielmoni, Antoine Houtain, Alessandro Derzelle, Karine Van Doninck, Jean-François Flot,
  BMC Bioinformatics 22:303, doi: https://doi.org/10.1186/s12859-021-04118-3


#### Requirements:

  - `sambamba`
  - `scipy`
  - `pandas`
  - `numpy`
  - `matplotlib`

#### Installation:

1. Install requirements:

```
# Install dependencies in a virtual environment (conda or virtualenv)
$ conda create -n happy-env pip numpy pandas scipy matplotlib docopt sambamba
$ conda activate happy-env
```

2. Install `HapPy`
- Using `pip` (previous version)
```
$ pip install happy-AntoineHo==0.2.1rc0
$ happy --help
```
- Using `git` (current version)
```
$ git clone https://github.com/AntoineHo/HapPy.git
$ python /path/to/happy/Hap.py --help
```

## 2. Main module
#### Usage:

```
$ python Hap.py -h
usage: Hap.py [-h] {coverage,estimate} ...

Estimate assembly haploidy based on base depth of coverage histogram.

positional arguments:
  {coverage,estimate}
    coverage            Compute coverage histogram for mapping file.
    estimate            Compute haploidy from coverage histogram.

optional arguments:
  -h, --help            show this help message and exit
```

## 3. Module coverage
This module runs `sambamba` on a read alignment file then reads the output depth file to obtain a coverage histogram.

#### Usage:
```
$ python Hap.py coverage -h
usage: Hap.py coverage [-h] [-t THREADS] [-d OUTDIR] MAP

positional arguments:
  MAP            <FILE> Sorted BAM file after mapping reads to the assembly.

optional arguments:
  -h, --help     show this help message and exit
  -t, --threads  <INT> Number of parallel threads allocated for sambamba. Default: 1
  -d, --outdir   <DIR> Path where the .cov and .hist files are written. Default: 'out'
```

## 4. Module estimate
Takes the .hist output file of module `coverage` and outputs metrics in a text file and optionnally as a graph. The size is provided with a value and a unit, ex: G for Gigabases, M for Megabases.

#### Usage:
```
$ python Hap.py estimate -h
usage: Hap.py estimate [-h] -s SIZE -o OUTSTATS [-m MIN_PEAK] [-p PROMINENCE]
                       [-vh VALLEY_HEIGHT] [-vp VALLEY_PROMINENCE] [-w WINDOW]
                       [-ll LIMIT_LOW] [-ld LIMIT_DIPLOID] [-lh LIMIT_HIGH]
                       [--plot] [--debug] [--no-smooth]
                       COV

positional arguments:
  COV                   <FILE> Coverage histogram file.

optional arguments:
  -h, --help               show this help message and exit
  -s, --size               <STRING> Estimated haploid genome size. (Recognized modifiers: K,M,G).
  -o, --outstats           <FILE> Path where haploidy value is written.
  -m, --min-peak           <INT> Minimum peak height (see SciPy doc). Default: 15000
  -p, --prominence         <INT> Minimum peak prominence (see SciPy docs). Default: 10000
  -vh, --valley-height     <INT> Minimum valley height (abs). Default: 15000
  -vp, --valley-prominence <INT> Minimum valley prominence (abs). Default: 10000
  -w, --window             <INT> Window length to use in savgol filter. Default: 41
  -ll, --limit-low         <INT> Lower threshold of coverage (lower than diploid peak). Default: auto.
  -ld, --limit-diploid     <INT> Middle threshold of coverage (between diploid and haploid peaks). Default: auto.
  -lh, --limit-high        <INT> Upper threshold of coverage (higher than haploid peak). Default: auto.
  --plot                   Generate histogram plot. Default: False
  --debug                  Generate more informative plots to debug. Default: False
  --no-smooth              Skip the smoothing step. Default: False
```

## 5. Example

Here is an example on how to use `HapPy`. `HapPy` requires a sorted BAM file as input. Here the PacBio long reads are mapped to the assembly with `minimap2`, and the output is sorted with `samtools`. The sorted BAM file is also indexed with `samtools`. The module depth computes the coverage histogram from the BAM file, and the module then estimates the haploidy metrics *H*. Here the max *x* value for the contaminant peak is set to 35, the max *x* value for the diploid peak is set to 120, and the size is set to 102 Mb.

```
# Align pacbio long reads on the assembly
$ minimap2 -ax map-pb assembly.fasta.gz pacbio_reads.fasta.gz --secondary=no | \
  samtools sort -o mapping_LR.map-pb.bam -T tmp.ali

# Index output BAM file
$ samtools index mapping_LR.map-pb.bam

# Obtain coverage histogram with sambamba
$ happy coverage -d happy_output mapping_LR.map-pb.bam

# Estimate Haploidy (manually input peak limits)
$ happy estimate --limit-low 35 --limit-diploid 120 --limit-high 200 -S 102M \
  -o happy_stats --plot happy_output/mapping_LR.map-pb.bam.hist
# Estimate Haploidy (try to detect peaks and limits automatically)
$ happy estimate -S 102M -p happy_stats --plot happy_output/mapping_LR.map-pb.bam.hist
```
