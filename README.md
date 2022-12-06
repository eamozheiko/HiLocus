# HiLocus - Translocation detection from capture Hi-C

## About HiLocus

This bioinformatic tool which help to find translocated regions from Hi-C data.
Key features of it algorims is possibility to detect translocation from sparse Hi-C matrix. For example: exome capture Hi-C.

HiLocus can detect both intra- and inter-chromosomal traslocations. It also finds the probability of translocation depending on the amount of region coverage, which helpful if coverage of you data is very heterogeneous.

HiLocus consist of two modules:
```
hilocus stat

hilocus trans

```

### hilocus stat

First module - ```hilocus stat``` obtain some key features of the quality of Hi-C libraries which may affect the accuracy of translocation detection.

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

Input for this module are files in .hic (highly recommended) or .pairs format:

```
    chromosome1 position1 chromosome2 position2
```

EXAMPLE: ```hilocus trans --case full/path/to/case --control full/path/to/control --binsize 10000 --name sample1 --outdir /path/to/outdir```

## Quick start

K562 and GM12878 (originally taken from GSE63525) are downsampled .hic files which can be used to test HiLocus trans. These .hic files only consist of chr3 - chr3 and chr3 - chr10 interactions and are only 1MB in resolution. You can find these files in the QuickStart folder. Use k562 as case, GM12878 as control and set binsize to 1000000. After running HiLocus trans, you can detect a FISH-validated reciprocal translocation with breakpoints in (chr3:48147000-48186000, chr10:86065000-86089000)

## Installation

### Dependencies

- Python >=3.9.15
- R >=4.2 (r-matrix, r-remotes, strawr (https://github.com/aidenlab/straw/tree/master/R))
- samtools >=1.6

### Install HiLocus

Install using conda

```
 conda install -c mozheiko -c bioconda hilocus
```

Manually installation 

1. Install HiLocus dependencies
2. Download HiLocus ```git clone https://github.com/eamozheiko/HiLocus.git```
3. Go to HiLocus directory, install it by ```python setup.py install```





