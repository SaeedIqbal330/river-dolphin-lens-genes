#!/bin/bash
# Step 2. Genome completeness with BUSCO (Table 1).
# Run in the BUSCO environment (environment_busco.yml).
H=${LENS_HOME:-$HOME/platanista_genomics}; cd "$H"

for key in 054084255 000442215 036417435 937001465 011762535 963924675; do
  g=$(find data/genomes -name "*${key}*.fna" | head -n 1)
  sp=$(echo "$g" | cut -d/ -f3)
  [ -d "busco_$sp" ] && { echo "$sp: done already"; continue; }
  busco -i "$g" -o "busco_$sp" -m genome -l cetartiodactyla_odb10 -c 8 --miniprot
done

grep -H "C:" busco_*/short_summary*.txt
