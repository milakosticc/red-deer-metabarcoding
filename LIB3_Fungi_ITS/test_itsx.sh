#!/usr/bin/env bash

set -euo pipefail

INPUT="test_dada2_out/test_ASVs.fasta"
OUTDIR="test_itsx_out"

mkdir -p "$OUTDIR"

ITSx \
  -i "$INPUT" \
  -o "$OUTDIR/fungi_test" \
  -t F \
  --save_regions ITS1 \
  --preserve T \
  --not_found T \
  --summary T \
  --positions T \
  --cpu 4

echo
echo "ITSx test completed."
echo "Output files:"
ls -lh "$OUTDIR"
