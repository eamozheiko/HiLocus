# HiLocus - Translocation detection from capture Hi-C

## About HiLocus

This bioinformatic tool which help to find translocated regions from Hi-C data.
Key features of it algorims is possibility to detect translocation from sparse Hi-C matrix. For example: exome capture Hi-C.

HiLocus can detect both intra- and inter-chromosomal traslocations. It also finds the probability of translocation depending on the amount of region coverage, which very helpful if coverage of you data is very heterogeneous

HiLocus consist of two modules:
```
hilocus stat

hilocus trans

```

### hilocus stat

First module - ```hilocus stat``` get some key features of quality of Hi-C libraries and prepare fastq files to ```hilocus trans``` and ```hilocus cnv``` input

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

EXAMPLE: ```hilocus stat --sam full/path/to/samfile --name sample1 --outdir full/path/to/outdir```

If your data in BAM or CRAM format, it is nessecary to convet it to SAM before runnig ```hilocus stat```

### hilocus trans

This is the main module which detect translocations from Hi-C data

Input for this module is filtred read pairs in the format:

```
    chromosome1 position1 chromosome2 position2
```

EXAMPLE: ```hilocus trans --case full/path/to/case_validpairs --control full/path/to/control_validpairs --binsize 20000 --name sample1 --outdir /path/to/outdir```


## Installation

### Dependencies

- Python 3.9
- R 3.6.1

### Install HiLocus
Manually installation 

1. Install HiLocus dependencies
2. Download HiLocus ```git clone https://github.com/eamozheiko/HiLocus.git```
3. Go to HiLocus directory, install it by ```python setup.py install```





