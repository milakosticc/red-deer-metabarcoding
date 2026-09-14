#!/bin/bash
set -euo pipefail

vsearch \
  --sintax taxonomy/otu_97/LIB3_OTU97_centroids.fasta \
  --db taxonomy/reference/UNITE_v10/utax_reference_dataset_19.02.2025.fasta \
  --tabbedout taxonomy/otu_97/LIB3_OTU97_UNITE_sintax_0.8.tsv \
  --sintax_cutoff 0.8 \
  --strand both \
  --threads 8 \
  --log taxonomy/otu_97/LIB3_OTU97_UNITE_sintax_0.8.log
