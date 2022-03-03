# ExoC - Translocation detection from Exome capture Hi-C

## About ExoC

This bioinformatic tool which help to find translocated regions from Hi-C data.
Key features of it algorims is possibility to detect translocation from sparce Hi-C matrix, for example, exome capture Hi-C.

ExoC can detect both intra- and inter-chromosomal traslocations. It also finds the probability of translocation depending on the amount of region coverage, which very helpful if coverage of you data is very ununiform

ExoC consist of three modules:

exoc stat

exoc trans

exoc cnv

## exoc stat

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

EXAMPLE: ```exoc stat -d /path/to/bamfile -b 20000 -n sample1 -o /path/to/outdir```

## exoc trans

This is the main module which detect translocations from Hi-C data

Input for this module is foltred read pairs in the format:

```
    chromosome1 position1 chromosome2 position2
```

EXAMPLE: ```exoc trans -i path/to/case_validpairs -c path/to/control_validpairs -b 20000 -n sample1 -o /path/to/outdir```

## exoc cnv

Such as detecting CNV from Hi-C data is has very low precision, so this module only calculate CNV pattern of known CNV regions

This can help the researcher to have more or less confidence in the CNV founded by standard methods 



