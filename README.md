# HapPy
**Hap**loidy with **Py**thon.

Easy haploidy and size completeness estimation.

## 1. General
This tool helps assess the haploidy and proximity to size completeness.

```
$ python HapPy/Hap.py -h
usage: Hap.py [-h] {depth,estimate} ...

Estimate assembly haploidy based on base depth of coverage histogram.

positional arguments:
  {depth,estimate}
    depth           Obtain a depth of coverage histogram for the assembly using sambamba depth base.
    estimate        Computes haploidy score based on the coverage distribution.

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
                       HIST OUT SIZE

positional arguments:
  HIST                  <STRING> A path to the histogram output of the `Hap.py depth` command (.hist file).
  OUT                   <STRING> A path for the output directory.
  SIZE                  <STRING> An expected assembly size (in bp) to compute the Total Size Score. Valid multipliers are (K, M, G) e.g.: 10K = 10000.
                                                                              
optional arguments:
  -h, --help            show this help message and exit
  -mc MAX_CONTAMINANT, --max-contaminant MAX_CONTAMINANT <INT> Maximum coverage of contaminants. Default: [35] 
  -md MAX_DIPLOID, --max-diploid MAX_DIPLOID <INT> Maximum coverage of the diploid peak. Default: [120]
  -mp MIN_PEAK, --min-peak MIN_PEAK <INT> Minimum peak height. Default: [150000]
  -p [], --plot []      Output plots. Default: False
```
