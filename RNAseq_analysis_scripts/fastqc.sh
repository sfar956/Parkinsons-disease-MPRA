#!/bin/bash -e

## Go to parent directory containing fastq files

cd ../data/rna_seq_sequencing_data/rnaseq_04122023/cleandata/

## perform FastQC analysis

for sample in *; 
do
    echo "Running fastqc for:" $sample;
    mkdir -p ../results/cleandata_fastqc/$sample;
    echo "Running fastqc on" $sample/*_R1.fq.gz
    fastqc -o ../results/cleandata_fastqc/$sample -t 20 -dir ../temp $sample/*_R1.fq.gz;
    echo "Running fastqc on" $sample/*_R2.fq.gz
    fastqc -o ../results/cleandata_fastqc/$sample -t 20 -dir ../temp $sample/*_R2.fq.gz;
    echo "Fastqc analysis completed for:" $sample "sample";
done

## compile FastQC results using multiQC

cd ../results/cleandata_fastqc

mkdir -p ../results/cleandata_fastqc/multiqc

multiqc -o ../results/cleandata_fastqc/multiqc .
