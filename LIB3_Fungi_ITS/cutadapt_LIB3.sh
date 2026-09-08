#!/bin/bash

set -euo pipefail

INPUT="working/Fungi_ITS/merged_raw"
OUTPUT="working/Fungi_ITS/primersCut_out"
LOGDIR="working/Fungi_ITS/logs"

mkdir -p "$OUTPUT"
mkdir -p "$LOGDIR"

FWD_PRIMER="GGAAGTAAAAGTCGTAACAAGG"
REV_PRIMER="GCTGCGTTCTTCATCGATGC"

for R1 in "$INPUT"/*_R1.fq.gz; do

    sample=$(basename "$R1" _R1.fq.gz)
    R2="$INPUT/${sample}_R2.fq.gz"

    echo "Processing $sample"

    if [ ! -f "$R2" ]; then
        echo "ERROR: Missing R2 for $sample"
        exit 1
    fi

    cutadapt \
        -g "^${FWD_PRIMER}" \
        -G "^${REV_PRIMER}" \
        --minimum-length 32 \
        --pair-filter=both \
        --cores=0 \
        -o "${OUTPUT}/${sample}_trimmed_R1.fastq.gz" \
        -p "${OUTPUT}/${sample}_trimmed_R2.fastq.gz" \
        "$R1" "$R2" \
        > "${LOGDIR}/${sample}_cutadapt.log" 2>&1

done
