# ExoC - Translocation detection from Exome capture Hi-C

## About ExoC

This bioinformatic tool which help to find translocated regions from Hi-C data.
Key features of it algorims is possibility to detect translocation from sparse Hi-C matrix. For example: exome capture Hi-C.

ExoC can detect both intra- and inter-chromosomal traslocations. It also finds the probability of translocation depending on the amount of region coverage, which very helpful if coverage of you data is very ununiform

ExoC consist of three modules:
```
exoc stat
```
```
exoc trans
```
```
exoc cnv
```

### exoc stat

First module - ```exoc stat``` get some key features of quality of Hi-C libraries and prepare fastq files to ```exoc trans``` and ```exoc cnv``` input

Key features of Hi-C quality:

- fracton of unique alignment in all pairs. (both reads alignes uniquely)
- fracton of multiple aligment in all pairs
- fracton of one side unmapped in all pairs
- fracton of both side unmapped in all pairs
- fracton of PCR duplicates in unique alignments
- fracton of dangling ends in unique alignments
- fracton of valid pairs in all pairs
- fracton of cis contacts in valid pairs
- fracton of trans contacts in valid pairs
- fracton of cis short range (distance less than 20kb) contacts in valid pairs
- fracton of cis long range (distance more than 20kb) contacts in valid pairs

EXAMPLE: ```exoc stat --sam /path/to/bamfile --name sample1 --outdir /path/to/outdir```

### exoc trans

This is the main module which detect translocations from Hi-C data

Input for this module is filtred read pairs in the format:

```
    chromosome1 position1 chromosome2 position2
```

EXAMPLE: ```exoc trans --case path/to/case_validpairs --control path/to/control_validpairs --binsize 20000 --name sample1 --outdir /path/to/outdir```

### exoc cnv

Such as detecting CNV from Hi-C data is has very low precision, so this module only calculate CNV pattern of known CNV regions

This can help the researcher to have more or less confidence in the CNV founded by standard methods 

EXAMPLE: ```exoc cnv --case path/to/case_validpairs --control path/to/control_validpairs --cnv /path/to/segmented_cnv.bed --name sample1 --outdir /path/to/outdir```

segmented CNV must be in .bed format:

```
    chromosome position_start position_end
```

## Installation

### Dependencies

- Python 3.9
- R 3.6.1

### Install ExoC
Manually installation 

1. Install ExoC dependencies
2. Download ExoC ```git clone https://github.com/eamozheiko/ExoC.git```
3. Go to ExoC directory, install it by ```python setup.py install```





