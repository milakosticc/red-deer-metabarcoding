#!/bin/bash
set -euo pipefail

vsearch \
--cluster_size taxonomy/otu_97/LIB2_ASVs_with_size.fasta \
--id 0.97 \
--sizein \
--sizeout \
--centroids taxonomy/otu_97/LIB2_OTU97_centroids_raw.fasta \
--uc taxonomy/otu_97/LIB2_OTU97_clusters.uc \
--threads 8 \
--log taxonomy/otu_97/LIB2_OTU97_vsearch.log
