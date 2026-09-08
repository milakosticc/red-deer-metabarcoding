#!/bin/bash

set -euo pipefail

INPUT="Cervus_elaphus/Plants_trnLc-trnLh"
OUTPUT="working/Plants_trnLc-trnLh/primersCut_out"
LOGDIR="working/Plants_trnLc-trnLh/logs"

mkdir -p "$OUTPUT"
mkdir -p "$LOGDIR"

FWD_PRIMER="CGAAATCGGTAGACGCTACG"
REV_PRIMER="CCATTGAGTCTCTGCACCTATC"

for d in "$INPUT"/LME*; do

    sample=$(basename "$d")

    R1=$(find "$d" -maxdepth 1 -type f -name "*_1.fq.gz")
    R2=$(find "$d" -maxdepth 1 -type f -name "*_2.fq.gz")

    echo "Processing $sample"

    cutadapt \
        -g "^${FWD_PRIMER}" \
        -G "^${REV_PRIMER}" \
        -o "${OUTPUT}/${sample}_trimmed_R1.fastq.gz" \
        -p "${OUTPUT}/${sample}_trimmed_R2.fastq.gz" \
        --minimum-length 32 \
        --pair-filter=both \
        --cores=0 \
        "$R1" "$R2" \
        > "${LOGDIR}/${sample}_cutadapt.log" 2>&1

done
