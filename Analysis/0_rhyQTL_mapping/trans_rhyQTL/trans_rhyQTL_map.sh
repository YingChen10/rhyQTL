#!/bin/bash

# all brain tissues
tissues=(
    "Brain-Amygdala"
    "Brain-Anteriorcingulatecortex_BA24"
    "Brain-Caudate_basalganglia"
    "Brain-CerebellarHemisphere"
    "Brain-Cerebellum"
    "Brain-Cortex"
    "Brain-FrontalCortex_BA9"
    "Brain-Hippocampus"
    "Brain-Hypothalamus"
    "Brain-Nucleusaccumbens_basalganglia"
    "Brain-Putamen_basalganglia"
    "Brain-Spinalcord_cervicalc-1"
    "Brain-Substantianigra"
)

# core clock gene
gene="NR1D1"


for tissue in "${tissues[@]}"; do
    echo "Processing $tissue with gene $gene..."
    Rscript 01_regerssion.R "$tissue" "$gene"
    Rscript 02_compare.R "$tissue" "$gene"
    Rscript 03_compare_multiple_times.R "$tissue" "$gene"
done

echo "All tissues completed."
