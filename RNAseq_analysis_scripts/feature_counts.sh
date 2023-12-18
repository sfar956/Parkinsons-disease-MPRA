#!/bin/bash -e

## Run featureCounts

featureCounts -a ../data/GRCh38.p13_reference_files/gencode.v43.primary_assembly.annotation.gtf \
    -o ../results/feature_counts/rnaseq_04122023_count_data.txt \
    ../results/star_mapping_outputs/*/*.sortedByCoord.out.bam \
    -F GTF \
    --largestOverlap \
    -p \
    --countReadPairs \
    -s 2 \
    -B \
    -C \
    -T 64 \
    -t gene \
    -g gene_name \
    --Rpath ../results/feature_counts/ \
    --tmpDir ../temp \
    --verbose
