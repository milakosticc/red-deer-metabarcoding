#!/bin/bash
set -euo pipefail

vsearch \
  --sintax taxonomy/otu_97/LIB2_OTU97_centroids.fasta \
  --db taxonomy/reference/trnL_CH_2026/Obitools_trnl_CH_SINTAX.fasta \
  --tabbedout taxonomy/otu_97/LIB2_OTU97_SINTAX_0.8.tsv \
  --sintax_cutoff 0.8 \
  --strand both \
  --threads 8 \
  --log taxonomy/otu_97/LIB2_OTU97_SINTAX_0.8.log
