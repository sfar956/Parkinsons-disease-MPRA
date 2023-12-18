#!/bin/bash -e

## reference genome indexing using STAR 

mkdir ../data/GRCh38.p13_indexed_genome

STAR --runThreadN 64 \
    --runMode genomeGenerate \
    --genomeDir ../data/GRCh38.p13_indexed_genome \
    --genomeFastaFiles ../data/GRCh38.p13_reference_files/GRCh38.primary_assembly.genome.fa \
    --sjdbGTFfile ../data/GRCh38.p13_reference_files/gencode.v43.primary_assembly.annotation.gtf \
    --sjdbOverhang 149

