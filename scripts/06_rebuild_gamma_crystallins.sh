#!/bin/bash
# Step 6. Rebuild CRYGC and CRYGB in P. gangetica from their own loci.
# The two genes lie side by side. In the genome-wide mapping, CRYGB's model ran into the CRYGC locus
# and CRYGC was left without a model. Each is rebuilt here by mapping its reference protein onto a
# window that holds only that gene. The rebuilt sequences are used in step 9 (selection tests)
# and step 11 (Table S2).
# Run after lens_pipeline.py map.
H=${LENS_HOME:-$HOME/platanista_genomics}; cd "$H/Lens_Project"
PLAT_DNA=$(find "$H/data/genomes" -name "*054084255*.fna" | head -n 1)

for job in "CRYGC JBLRPT010000008.1:30248000-30254000" "CRYGB JBLRPT010000008.1:30236000-30244000"; do
  set -- $job; G=$1; REGION=$2; g=$(echo "$G" | tr 'A-Z' 'a-z')
  samtools faidx "$PLAT_DNA" "$REGION" > ${g}_window.fa
  python3 - "$G" > ${g}_ref.faa << 'PYEOF'
import sys
seqs, n = {}, None
for l in open("refs.faa"):
    l = l.strip()
    if l.startswith(">"): n = l[1:]; seqs[n] = ""
    else: seqs[n] += l
print(f">{sys.argv[1]}\n{seqs[sys.argv[1]]}")
PYEOF
  miniprot --gff ${g}_window.fa ${g}_ref.faa 2>/dev/null > ${g}.gff
  grep -v "^##PAF" ${g}.gff > ${g}_clean.gff
  gffread ${g}_clean.gff -g ${g}_window.fa -y ${g}_prot.faa -x ${g}_cds.fa
  echo "== $G (window $REGION)"
  grep -P "\tmRNA\t" ${g}.gff | cut -f9 | head -n 1
  grep -P "\tCDS\t" ${g}.gff | awk '{print "   exon", $4"-"$5, "(window coordinates)"}'
  cat ${g}_prot.faa
done
