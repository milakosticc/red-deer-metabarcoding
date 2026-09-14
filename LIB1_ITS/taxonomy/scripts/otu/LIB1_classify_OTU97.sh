#!/bin/bash
set -euo pipefail

vsearch \
--sintax taxonomy/otu_97/LIB1_OTU97_centroids.fasta \
--db taxonomy/reference/BOLD_plant_ITS/BOLD_Plantae_ITS2_SINTAX.fasta \
--tabbedout taxonomy/otu_97/LIB1_OTU97_BOLD_sintax_0.8.tsv \
--sintax_cutoff 0.8 \
--strand both \
--threads 8 \
--log taxonomy/otu_97/LIB1_OTU97_BOLD_sintax_0.8.log
