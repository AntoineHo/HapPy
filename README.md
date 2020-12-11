# HapPy
**Hap**loidy using **Py**thon.

Easy haploidy estimation.

[![DOI](https://zenodo.org/badge/299235590.svg)](https://zenodo.org/badge/latestdoi/299235590)

## 1. General
This tool assesses the haploidy *H* of a given assembly.
*H* is defined as the fraction of the bases of the genome that are in the collapsed peak *C*. This metrics is calculated as *H*=*C*/(*C*+*U*/2), where *C* is the size (area) of the collapsed peak and *U* the size of the uncollapsed peak in the per-base coverage histogram of the assembly.

For more information, see:
  **Overcoming uncollapsed haplotypes in long-read assemblies of non-model organisms**, 
  Nadège Guiglielmoni, Antoine Houtain, Alessandro Derzelle, Karine Van Doninck, Jean-François Flot,
  bioRxiv 2020, doi: https://doi.org/10.1101/2020.03.16.993428 

### Requirements: 

- `sambamba`
- `scipy`
- `pandas`
- `numpy`

```
$ python HapPy/main.py -h
Estimate assembly haploidy based on base depth of coverage histogram.

usage:
    HapPy [-hv] <command> [<args>...]

options:
    -h, --help                  shows the help
    -v, --version               shows the version

The subcommands are:
    coverage    Compute coverage histogram.
    estimate    Finds peaks and modality, then computes scores of haploidy.
```

## 2. Module coverage
This module runs `sambamba` on a read alignment file then reads the output depth file to obtain a coverage histogram.

```
$ python HapPy/main.py coverage -h

Coverage histogram command
    Compute coverage histogram for mapping file.

    usage:
        coverage [--threads=1] --outdir=DIR <mapping.bam>
        
    arguments:
        mapping.bam              Sorted BAM file after mapping reads to the assembly.
        
    options:
        -t, --threads=INT        Number of parallel threads allocated for 
                                 sambamba [default: 1].
        -d, --outdir=DIR         Path where the .cov and .hist files are written.
```

## 3. Module estimate
Takes the .hist output file of module `coverage` and outputs metrics in a text file and optionnally as a graph. The size is provided with a value and a unit, ex: G for Gigabases, M for Megabases.

### Usage:
```
$ python HapPy/main.py estimate -h 
Estimate command
    Compute haploidy from coverage histogram.

    usage:
        estimate [--max-contaminant=INT] [--max-diploid=INT] [--min-peak=INT] 
                 --size=INT --outstats=FILE [--plot] <coverage.hist>
        
    arguments:
        coverage.hist               Coverage histogram.
        
    options:
        -C, --max-contaminant=INT   Maximum coverage of contaminants.
        -D, --max-diploid=INT       Maximum coverage of the diploid peak.
        -M, --min-peak=INT          Minimum peak height.
        -S, --size=INT              Estimated haploid genome size.
        -O, --outstats=FILE         Path where the AUC ratio and TSS values are written.
        -p, --plot                  Generate histogram plot.
```

## 4. Example

Here is an example on how to use `HapPy`. `HapPy` requires a sorted BAM file as input. Here the PacBio long reads are mapped to the assembly with `minimap2`, and the output is sorted with `samtools`. The sorted BAM file is also indexed with `samtools`. The module depth computes the coverage histogram from the BAM file, and the module then estimates the haploidy metrics *H*. Here the max *x* value for the contaminant peak is set to 35, the max *x* value for the diploid peak is set to 120, and the size is set to 102 Mb.

```
minimap2 -ax map-pb assembly.fasta.gz pacbio_reads.fasta.gz --secondary=no | samtools sort -o mapping_LR.map-pb.bam -T tmp.ali
samtools index mapping_LR.map-pb.bam
main.py coverage -d happy_output mapping_LR.map-pb.bam 
Hap.py estimate --max-contaminant 35 --max-diploid 120 -S 102M -O happy_stats.txt -p happy_output/mapping_LR.map-pb.bam.hist
```
