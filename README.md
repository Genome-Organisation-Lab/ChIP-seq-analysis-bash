# ChIP-seq analysis

A bash script for paired-end ChIP-seq data processing and peak-calling, based on the established [REMAP pipeline](https://doi.org/10.1093/nar/gkab996).

## Pre-requisites and dependencies
* [python3-biopython](https://github.com/biopython/biopython)
* [sra-toolkit](https://github.com/ncbi/sra-tools)
* [trim-galore](https://github.com/FelixKrueger/TrimGalore)
* [bowtie2](https://github.com/BenLangmead/bowtie2)
* [samtools](https://github.com/samtools/samtools)
* [bedtools](https://github.com/arq5x/bedtools2)
* [MACS2](https://github.com/macs3-project/MACS)

## Usage
```
$ bash run_analysis.sh {Test_Rep1} {Control_Rep1} {Test_Rep2} {Control_Rep2}
```
Inputs should not include .fastq extension. 

Inclusion of multiple replicates and controls should be specified in the script under

> ##Read inputs


## Customisation

The pipeline runs in 5 main stages, each of which needs to be configured based on the input specifications.

## Output
For each input, MACS2 putputs two directories:
* {Test_Rep}peakCallingResultsNarrow contains narrow peaks files
* {Test_Rep}peakCallingResultsBroad contains broad peaks files

Outputs from multiple replicates will be merged.
