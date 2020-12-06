# HapPy
**Hap**loidy using **Py**thon.

Easy haploidy estimation.

[![DOI](https://zenodo.org/badge/299235590.svg)](https://zenodo.org/badge/latestdoi/299235590)

## 1. General
This tool helps assess the haploidy H of a given assembly.
H is defined as the fraction of the bases of the genome that are in the collapsed peak C. This metrics is calculated as H=C/(C+U/2), where C is the size of the collapsed peak and U the size of the uncollapsed peak in the per-base coverage histogram of the assembly.

For more information, see:
  **Overcoming uncollapsed haplotypes in long-read assemblies of non-model organisms**, 
  Nadège Guiglielmoni, Antoine Houtain,Alessandro Derzelle, Karine van Doninck, Jean-François Flot,
  bioRxiv 2020, doi: https://doi.org/10.1101/2020.03.16.993428 

### Requirements: 

- sambamba
- scipy
- pandas
- numpy

```
$ python HapPy/Hap.py -h
usage: Hap.py [-h] {depth,estimate} ...

Estimate assembly haploidy based on coverage depth.

positional arguments:
  {depth,estimate}
    depth           Obtain a depth of coverage histogram for the assembly using sambamba depth base.
    estimate        Computes the haploidy score based on the coverage distribution.

optional arguments:
  -h, --help        show this help message and exit
```

## 2. Module depth
This module runs sambamba on a read alignment file then reads the output depth file to obtain a coverage histogram.

```
$ python HapPy/Hap.py depth -h
usage: Hap.py depth [-h] [-t THREADS] BAM OUT

positional arguments:
  BAM                   <STRING> A path to the reads alignment on the assembly (.bam file).
  OUT                   <STRING> An output path for the coverage files (.cov and .hist files).                                     

optional arguments:
  -h, --help            show this help message and exit
  -t THREADS, --threads THREADS <INT> Number of threads for sambamba. Default: [4].
```

## 3. Module estimate
Takes .hist output file of `Hap.py depth` and outputs metrics in a text file and optionnally a graph.

### Usage:
```
$ python HapPy/Hap.py estimate -h 
usage: Hap.py estimate [-h] [-mc MAX_CONTAMINANT] [-md MAX_DIPLOID]
                       [-mp MIN_PEAK] [-p ]
                       HIST OUT

positional arguments:
  HIST                  <STRING> A path to the histogram output of the `Hap.py depth` command (.hist file).
  OUT                   <STRING> A path for the output directory.
                                                                              
optional arguments:
  -h, --help            show this help message and exit
  -mc MAX_CONTAMINANT, --max-contaminant MAX_CONTAMINANT <INT> Maximum coverage of contaminants. Default: [35] 
  -md MAX_DIPLOID, --max-diploid MAX_DIPLOID <INT> Maximum coverage of the diploid peak. Default: [120]
  -mp MIN_PEAK, --min-peak MIN_PEAK <INT> Minimum peak height. Default: [150000]
  -p, --plot      Output plots. Default: False
```

## 4. Example

Here is an example on how to use HapPy. HapPy requires a sorted BAM file as input. Here the PacBio long reads are mapped to the assembly with minimap2, and the output is sorted with samtools. The sorted BAM file is also indexed with samtools. The module depth computes the coverage histogram from the BAM file, and then the module estimate computes the haploidy metrics H. Here the max x value for the contaminant peak is set to 35, the max x value for the diploid peak is set to 120, and the min y for a peak is set to 150000. The expected genome size is set to 102M. 

```
minimap2 -ax map-pb assembly.fasta.gz pacbio_reads.fasta.gz --secondary=no | samtools sort -o mapping_LR.map-pb.bam -T tmp.ali
samtools index mapping_LR.map-pb.bam
Hap.py depth mapping_LR.map-pb.bam happy_output
Hap.py estimate --max-contaminant 35 --max-diploid 120 --min-peak 150000 happy_output/mapping_LR.map-pb.bam.hist . 
```
