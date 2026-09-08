#!/bin/bash

set -euo pipefail

INPUT="Cervus_elaphus/Fungi_ITS"
OUTPUT="working/Fungi_ITS/merged_raw"

mkdir -p "$OUTPUT"

for dir in "$INPUT"/*; do

    sample=$(basename "$dir")

    mapfile -t R1_FILES < <(
        find "$dir" -maxdepth 1 -type f -name "*_1.fq.gz" | sort
    )

    mapfile -t R2_FILES < <(
        find "$dir" -maxdepth 1 -type f -name "*_2.fq.gz" | sort
    )

    echo "$sample: R1=${#R1_FILES[@]} R2=${#R2_FILES[@]}"

    if [ "${#R1_FILES[@]}" -ne "${#R2_FILES[@]}" ]; then
        echo "ERROR: unequal number of R1 and R2 files for $sample"
        exit 1
    fi

    if [ "${#R1_FILES[@]}" -eq 1 ]; then

        cp "${R1_FILES[0]}" \
           "$OUTPUT/${sample}_R1.fq.gz"

        cp "${R2_FILES[0]}" \
           "$OUTPUT/${sample}_R2.fq.gz"

    elif [ "${#R1_FILES[@]}" -gt 1 ]; then

        zcat "${R1_FILES[@]}" | gzip -c \
            > "$OUTPUT/${sample}_R1.fq.gz"

        zcat "${R2_FILES[@]}" | gzip -c \
            > "$OUTPUT/${sample}_R2.fq.gz"

    else

        echo "ERROR: no FASTQ files found for $sample"
        exit 1

    fi

done
