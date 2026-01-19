# ChIP-seq analysis

A bash script for paired-end ChIP-seq data processing and peak-calling, based on the established [REMAP pipeline](https://doi.org/10.1093/nar/gkab996).

## Pre-requisites and dependencies
* [biopython 1.81](https://github.com/biopython/biopython)
* [sra-toolkit 2.10.9](https://github.com/ncbi/sra-tools)
* [trim-galore 8.7.32](https://github.com/FelixKrueger/TrimGalore)
* [bowtie2 2.5.1](https://github.com/BenLangmead/bowtie2)
* [samtools 1.18](https://github.com/samtools/samtools)
* [bedtools 2.31.0](https://github.com/arq5x/bedtools2)
* [MACS2 2.2.9.1](https://github.com/macs3-project/MACS)

## Usage
```
$ bash run_analysis.sh Test_Rep1 ( Control_Rep1 Test_Rep2 Control_Rep2 )
```
Inputs should not include .fastq extension. 

## Customisation

The pipeline runs in 5 main stages, each of which needs to be configured based on the input specifications.

*fastaLength* is calculated using
```
$ getFastaLength.py genomeFasta
```

## Output
For each input, MACS2 outputs two directories:
* {Test_RepN}peakCallingResultsNarrow contains narrow peaks files
* {Test_RepN}peakCallingResultsBroad contains broad peaks files

Outputs from multiple replicates will be merged.
