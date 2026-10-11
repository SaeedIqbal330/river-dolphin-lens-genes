#!/bin/bash
# Step 1. Download the six genome assemblies from NCBI.
# Needs the NCBI datasets tool (conda install -c conda-forge ncbi-datasets-cli).
# Genomes go to $LENS_HOME/data/genomes/<Species>/ (LENS_HOME defaults to ~/platanista_genomics).
set -e
H=${LENS_HOME:-$HOME/platanista_genomics}

declare -A GENOMES=(
  [Platanista_gangetica]=GCA_054084255.1
  [Lipotes_vexillifer]=GCF_000442215.2
  [Inia_geoffrensis]=GCA_036417435.1
  [Orcinus_orca]=GCF_937001465.1
  [Tursiops_truncatus]=GCA_011762535.1
  [Phocoena_phocoena]=GCF_963924675.1
)

for sp in "${!GENOMES[@]}"; do
  acc=${GENOMES[$sp]}
  mkdir -p "$H/data/genomes/$sp" && cd "$H/data/genomes/$sp"
  if ls ncbi_dataset/data/$acc/*.fna > /dev/null 2>&1; then
    echo "$sp: already downloaded"; continue
  fi
  datasets download genome accession "$acc" --include genome --filename genome.zip
  unzip -o -q genome.zip && rm genome.zip
  f=$(ls ncbi_dataset/data/$acc/*.fna)
  echo "$sp: $f ($(du -h "$f" | cut -f1), $(grep -c '>' "$f") sequences)"
done
