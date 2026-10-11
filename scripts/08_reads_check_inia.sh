#!/bin/bash
# Step 8. Check the CRYBA1 frameshift in I. geoffrensis in raw reads (Table 3, Table S5).
#   SRR11430609 = Illumina, separate individual (60 million reads streamed)
#   SRR30296906 = PacBio HiFi, assembly project (400,000 reads streamed)
# An indel is counted only if it appears in at least 30% of the reads at a position, so that
# occasional sequencing errors in single reads do not count. With deep coverage, a heterozygous
# frameshift would show an indel in about half the reads.
# Takes a few hours, mostly download. Run in the background, e.g.:
#   nohup bash 08_reads_check_inia.sh > reads_inia.txt 2>&1 &
H=${LENS_HOME:-$HOME/platanista_genomics}
mkdir -p "$H/Lens_reads/Inia" && cd "$H/Lens_reads/Inia"
INIA_DNA=$(find "$H/data/genomes" -name "*036417435*.fna" | head -n 1)

# 100-kb window around CRYBA1; frameshift at window position 50,001, first stop at 48,892
if [ ! -s window.fa ]; then
  samtools faidx "$INIA_DNA" CM070917.1:21079508-21179508 | sed 's/^>.*/>CRYBA1/' > window.fa
  samtools faidx window.fa
fi

# run, reads to stream, minimap2 preset
for job in "SRR11430609 60000000 sr" "SRR30296906 400000 map-hifi"; do
  set -- $job; RUN=$1; N=$2; PRESET=$3
  if [ ! -s ${RUN}.bam ]; then
    fastq-dump -X $N --stdout $RUN \
      | minimap2 -ax $PRESET window.fa - 2>/dev/null \
      | samtools view -b -F 4 -q 20 - | samtools sort -o ${RUN}.bam
    samtools index ${RUN}.bam
  fi
done

python3 - << 'PYEOF'
import subprocess, re
sites = [("CRYBA1", 50001, "frameshift region"), ("CRYBA1", 48892, "first stop codon")]
labels = {"SRR11430609": "Illumina, separate individual", "SRR30296906": "PacBio, assembly project"}
for run, label in labels.items():
    print(f"\n=== {run}  ({label})")
    for win, pos, what in sites:
        out = subprocess.run(["samtools", "mpileup", "-f", "window.fa", "-r", f"{win}:{pos-60}-{pos+60}", f"{run}.bam"],
                             capture_output=True, text=True).stdout.splitlines()
        depths, indels, mism = [], [], []
        for l in out:
            f = l.split("\t"); p, d, bases = int(f[1]), int(f[3]), f[4]; depths.append(d)
            n_ind, n_mis, i = 0, 0, 0
            while i < len(bases):
                c = bases[i]
                if c == "^": i += 2; continue
                if c in "+-":
                    m = re.match(r"[+-](\d+)", bases[i:]); n_ind += 1
                    i += len(m.group(0)) + int(m.group(1)); continue
                if c in "ACGTNacgtn": n_mis += 1
                i += 1
            if n_ind and d and n_ind / d >= 0.3: indels.append(f"{p}({n_ind}/{d})")
            if n_mis and d and n_mis / d >= 0.5: mism.append(f"{p}({n_mis}/{d})")
        cov = f"depth {min(depths)}-{max(depths)}" if depths else "NO READS"
        verdict = "reads AGREE with the assembly" if depths and min(depths) >= 2 and not indels and not mism else "CHECK"
        print(f"  {what:18} {cov:14}  indels in >=30% of reads: {', '.join(indels) or 'none'}   "
              f"majority mismatches: {', '.join(mism) or 'none'}   -> {verdict}")
PYEOF
