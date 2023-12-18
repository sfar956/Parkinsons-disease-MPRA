#!/bin/bash -e

## Read mapping using STAR

cd ../data/rna_seq_sequencing_data/rnaseq_04122023/cleandata/

for dir in *;
do
    echo $dir;
    STAR --runMode alignReads \
        --runThreadN 64 \
        --genomeDir ../data/GRCh38.p13_indexed_genome \
        --readFilesIn $dir/*_R1.fq.gz $dir/*_R2.fq.gz \
        --readFilesCommand zcat \
        --sjdbGTFfile ../data/GRCh38.p13_reference_files/gencode.v43.primary_assembly.annotation.gtf \
        --sjdbOverhang 149 \
        --runRNGseed 444 \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode TranscriptomeSAM \
        --outFileNamePrefix ../results/star_mapping_outputs/$dir/$dir \
        --outTmpDir ../temp;
done
