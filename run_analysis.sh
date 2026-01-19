#!/bin/bash

##Processes ChIP-seq data (paired-end) and performs peak-calling with macs2, approximating the established REMAP pipeline (https://doi.org/10.1093/nar/gkab996).
##For the nextflow version of the analysis pipeline, check the repository "ChIP-seq-analysis-nextflow" (https://github.com/Genome-Organisation-Lab/ChIP-seq-analysis-nextflow).
#Usage: bash run_analysis.sh $1 $2

##Pre-requisites and dependencies
##sra-toolkit trim-galore bowtie2 samtools bedtools python3-biopython macs2 

########################################################################################

##Read sequence inputs
##comment out unused inputs
Sample_Rep1=$1
Control_Rep1=$2
Sample_Rep2=$3
Control_Rep2=$4

inDir=/path/to/input/directory
outDir=/path/to/output/directory
mkdir -p $outDir

##Set inputs and parameters

trimCORES=1 ##change based on core availability
CORES=1 ##change based on core availability

#Reference genome and basic information
genome=GCA_900044135.1_GZPH1RResV1_genomic.fna
index=PH1
#fastaLength=$(getFastaLength.py !{genomeFasta})
fastaLength=59900000 ##analysis from the script
epiMARK=H3K27me3 ##epigenetic mark under study
genotype=fusarium ##genotype background under study

########################################################################################

## Run pipeline
##0## Create epigenetic mark specific library
mkdir -p $genotype_$epiMARK

########################################################################################

##1## Building bowtie2 index for reference genome

bowtie2-build --threads $CORES $refDir/$genome $index

########################################################################################

##2## Trim the raw "Sample" and "Control" ChIP-seq reads using trim galore: paired end

trim_galore --output_dir $outDir/ --length 30 --quality 20 --stringency 1 -e 0.1 --paired -j $trimCORES $inDir/${Sample_Rep1}_1.fastq $inDir/${Sample_Rep1}_2.fastq

## Un-comment the following commands if you have "control"
#trim_galore --output_dir $outDir/ --length 30 --quality 20 --stringency 1 -e 0.1 --paired -j $trimCORES $inDir/${Control_Rep1}_1.fastq $inDir/${Control_Rep1}_2.fastq

## Un-comment the following commands if you have "replicates"
#trim_galore --output_dir $outDir/ --length 30 --quality 20 --stringency 1 -e 0.1 --paired -j $CORES $inDir/${Sample_Rep2}_1.fastq $inDir/${Sample_Rep2}_2.fastq
#trim_galore --output_dir $outDir/ --length 30 --quality 20 --stringency 1 -e 0.1 --paired -j $CORES $inDir/${Control_Rep2}_1.fastq $inDir/${Control_Rep2}_2.fastq

echo "Trimming finished for all samples"

########################################################################################

##3## Align the QC test and control reads against the genome using bowtie2

bowtie2 --threads $CORES --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 X 700 -x $refDir/$index -1 $outDir/${Sample_Rep1}_1_val_1.fq -2 $outDir/${Sample_Rep1}_2_val_2.fq -S $outDir/$Sample_Rep1.sam 
samtools view --threads $CORES -bhS $outDir/${Sample_Rep1}.sam -o $outDir/${Sample_Rep1}.bam
samtools index --threads $CORES $outDir/${Sample_Rep1}.bam

## Un-comment the following commands if you have "control"
#bowtie2 --threads $CORES --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 X 700 -x $refDir/$index -1 $outDir/${Control_Rep1}_1_val_1.fq -2 $outDir/${Control_Rep1}_2_val_2.fq -S $outDir/$Control_Rep1.sam 
#samtools view --threads $CORES -bhS $outDir/${Control_Rep1}.sam -o $outDir/${Control_Rep1}.bam
#samtools index --threads $CORES $outDir/${Control_Rep1}.bam

## Un-comment the following commands if you have "replicates"
#bowtie2 --threads $CORES --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 X 700 -x $refDir/$index -1 $outDir/${Sample_Rep2}_1_val_1.fq -2 $outDir/${Sample_Rep2}_2_val_2.fq -S $outDir/$Sample_Rep2.sam 
#samtools view --threads $CORES -bhS $outDir/${Sample_Rep2}.sam -o $outDir/${Sample_Rep2}.bam
#samtools index --threads $CORES $outDir/${Sample_Rep2}.bam

#bowtie2 --threads $CORES --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 X 700 -x $refDir/$index -1 $outDir/${Control_Rep2}_1_val_1.fq -2 $outDir/${Control_Rep2}_2_val_2.fq -S $outDir/$Control_Rep2.sam 
#samtools view --threads $CORES -bhS $outDir/${Control_Rep2}.sam -o $outDir/${Control_Rep2}.bam
#samtools index --threads $CORES $outDir/${Control_Rep2}.bam

echo "Mapping finished for all samples"

########################################################################################

##4## Use samtools to remove PCR duplicates

samtools sort --threads $CORES -n $outDir/${Sample_Rep1}.sam -o $outDir/${Sample_Rep1}_sorted.sam 
samtools fixmate --threads $CORES -m $outDir/${Sample_Rep1}_sorted.sam $outDir/${Sample_Rep1}_sorted.scored.sam  # m - add mate tag
samtools sort --threads $CORES $outDir/${Sample_Rep1}_sorted.scored.sam -o $outDir/${Sample_Rep1}_sorted.scored.sorted.sam
samtools markdup --threads $CORES -O BAM -r $outDir/${Sample_Rep1}_sorted.scored.sorted.sam $outDir/${Sample_Rep1}_deduplicated.bam #r - remove duplicate reads
samtools index --threads $CORES $outDir/${Sample_Rep1}_deduplicated.bam

## Un-comment the following commands if you have "control"
#samtools sort --threads $CORES -n $outDir/${Control_Rep1}.sam -o $outDir/${Control_Rep1}_sorted.sam 
#samtools fixmate --threads $CORES -m $outDir/${Control_Rep1}_sorted.sam $outDir/${Control_Rep1}_sorted.scored.sam 
#samtools sort --threads $CORES $outDir/${Control_Rep1}_sorted.scored.sam -o $outDir/${Control_Rep1}_sorted.scored.sorted.sam
#samtools markdup --threads $CORES -O BAM -r $outDir/${Control_Rep1}_sorted.scored.sorted.sam $outDir/${Control_Rep1}_deduplicated.bam
#samtools index --threads $CORES $outDir/${Sample_Rep2}_deduplicated.bam

## Un-comment the following commands if you have "replicates"
#samtools sort --threads $CORES -n $outDir/${Sample_Rep2}.sam -o $outDir/${Sample_Rep2}_sorted.sam 
#samtools fixmate --threads $CORES -m $outDir/${Sample_Rep2}_sorted.sam $outDir/${Sample_Rep2}_sorted.scored.sam 
#samtools sort --threads $CORES $outDir/${Sample_Rep2}_sorted.scored.sam -o $outDir/${Sample_Rep2}_sorted.scored.sorted.sam
#samtools markdup --threads $CORES -O BAM -r $outDir/${Sample_Rep2}_sorted.scored.sorted.sam $outDir/${Sample_Rep2}_deduplicated.bam
#samtools index --threads $CORES $outDir/${Control_Rep1}_deduplicated.bam

#samtools sort --threads $CORES -n $outDir/${Control_Rep2}.sam -o $outDir/${Control_Rep2}_sorted.sam 
#samtools fixmate --threads $CORES -m $outDir/${Control_Rep2}_sorted.sam $outDir/${Control_Rep2}_sorted.scored.sam 
#samtools sort --threads $CORES $outDir/${Control_Rep2}_sorted.scored.sam -o $outDir/${Control_Rep2}_sorted.scored.sorted.sam
#samtools markdup --threads $CORES -O BAM -r $outDir/${Control_Rep2}_sorted.scored.sorted.sam $outDir/${Control_Rep2}_deduplicated.bam
#samtools index --threads $CORES $outDir/${Control_Rep2}_deduplicated.bam

echo "PCR duplicates removed"

########################################################################################

##5## Call both broad and narrow peaks using the MACS2 peakcaller

macs2 callpeak --broad -g $fastaLength -q 0.00001 -t $outDir/${Sample_Rep1}_deduplicated.bam --outdir $outDir --name ${Sample_Rep1}_broad -f BAMPE  --broad-cutoff 0.00001
macs2 callpeak -g $fastaLength -q 0.00001 -t $outDir/${Sample_Rep1}_deduplicated.bam --outdir $outDir --name ${Sample_Rep1}_narrow -f BAMPE

## Un-comment the following commands if you have "control"
#macs2 callpeak -c $outDir/${Control_Rep1}_deduplicated.bam --broad -g $fastaLength -q 0.00001 -t $outDir/${Sample_Rep1}_deduplicated.bam --outdir $outDir --name ${Sample_Rep1}_broad -f BAMPE  --broad-cutoff 0.00001
#macs2 callpeak -c $outDir/${Control_Rep1}_deduplicated.bam -g $fastaLength -q 0.00001 -t $outDir/${Sample_Rep1}_deduplicated.bam --outdir $outDir --name ${Sample_Rep1}_broad -f BAMPE

## Un-comment the following commands if you have "replicates"
#macs2 callpeak --broad -g $fastaLength -q 0.00001 -t $outDir/${Sample_Rep2}_deduplicated.bam --outdir $outDir --name ${Sample_Rep2}_broad -f BAMPE  --broad-cutoff 0.00001
#macs2 callpeak -g $fastaLength -q 0.00001 -t $outDir/${Sample_Rep2}_deduplicated.bam --outdir $outDir --name ${Sample_Rep2}_narrow -f BAMPE

#macs2 callpeak -c $outDir/${Control_Rep2}_deduplicated.bam --broad -g $fastaLength -q 0.00001 -t $outDir/${Sample_Rep2}_deduplicated.bam --outdir $outDir --name ${Sample_Rep2}_broad -f BAMPE  --broad-cutoff 0.00001
#macs2 callpeak -c $outDir/${Control_Rep2}_deduplicated.bam -g $fastaLength -q 0.00001 -t $outDir/${Sample_Rep2}_deduplicated.bam --outdir $outDir --name ${Sample_Rep2}_broad -f BAMPE

echo "Peak calling finished"

## Merge replicates 

cat $outDir/*_peaks.narrowPeak > $outDir/unmergedNarrowPeaks.bed
sort -k1,1 -k2,2n $outDir/unmergedNarrowPeaks.bed > $outDir/unmergedNarrowPeaks.sorted.bed
bedtools merge -i $outDir/unmergedNarrowPeaks.sorted.bed > $outDir/${genotype}_${epiMARK}_mergedNarrowPeaks.bed

cat $outDir/*_peaks.broadPeak > $outDir/unmergedBroadPeaks.bed
sort -k1,1 -k2,2n $outDir/unmergedBroadPeaks.bed > $outDir/unmergedBroadPeaks.sorted.bed
bedtools merge -i $outDir/unmergedBroadPeaks.sorted.bed > $outDir/${genotype}_${epiMARK}_mergedBroadPeaks.bed
    
#rm $outDir/*.sam
#rm $outDir/*.fq
#rm $outDir/*_trimming_report.txt
    
echo "peak_calling pipeline finished for ${genotype}_${epiMARK}"
