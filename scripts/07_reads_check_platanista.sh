#!/bin/bash
# Step 7. Check every P. gangetica break in raw sequencing reads (Table 3, Table S5).
# Reads are streamed from the SRA (fastq-dump) and mapped with minimap2 onto 100-kb windows
# around each break; only reads with mapping quality >= 20 are kept.
#   SRR13005646 = second individual;  SRR39498552 = the assembly's own animal.
# A break is accepted when at least two reads cross it, no read shows an insertion or deletion
# there, and no majority of reads changes a stop codon.
# Takes a few hours, mostly download. Run in the background, e.g.:
#   nohup bash 07_reads_check_platanista.sh > reads_platanista.txt 2>&1 &
H=${LENS_HOME:-$HOME/platanista_genomics}
PLAT_DNA=$(find "$H/data/genomes" -name "*054084255*.fna" | head -n 1)
RUNS="SRR13005646 SRR39498552"

# ---------------------------------------------------------------- Part A: CRYGD
# 100-kb window JBLRPT010000008.1:30210000-30310000; the frameshifted exon sits at window 53,712-53,983.
mkdir -p "$H/CRYGD_reads" && cd "$H/CRYGD_reads"
if [ ! -s crygd_win.fa ]; then
  samtools faidx "$PLAT_DNA" JBLRPT010000008.1:30210000-30310000 | sed 's/^>.*/>crygd_win/' > crygd_win.fa
  samtools faidx crygd_win.fa
fi
for RUN in $RUNS; do
  if [ ! -s ${RUN}.bam ]; then
    fastq-dump -X 20000000 --stdout $RUN \
      | minimap2 -ax sr crygd_win.fa - 2>/dev/null \
      | samtools view -b -F 4 -q 20 - | samtools sort -o ${RUN}.bam
    samtools index ${RUN}.bam
  fi
  echo "== CRYGD, $RUN: mean read depth over the frameshifted exon"
  samtools depth -a -r crygd_win:53712-53983 ${RUN}.bam | awk '{s+=$3} END {print s/NR}'
  echo "== CRYGD, $RUN: differences between reads and the assembly in the coding exons"
  bcftools mpileup -f crygd_win.fa -r crygd_win:51200-54000 ${RUN}.bam 2>/dev/null \
    | bcftools call -mv 2>/dev/null | grep -v "^#" | cut -f2,4,5,10
done

# ---------------------------------------------------------------- Part B: CRYBB2, CRYGB, BFSP1, BFSP2
# 100-kb windows; each break sits at window position 50,001 (BFSP1's stop at 50,394).
mkdir -p "$H/Lens_reads" && cd "$H/Lens_reads"
regions="CRYBB2:JBLRPT010000015.1:73195827-73295827 CRYGB:JBLRPT010000008.1:30190533-30290533 BFSP1:JBLRPT010000016.1:51214043-51314043 BFSP2:JBLRPT010000001.1:15260985-15360985"
if [ ! -s windows.fa ]; then
  for r in $regions; do
    name=${r%%:*}; reg=${r#*:}
    samtools faidx "$PLAT_DNA" "$reg" | sed "s/^>.*/>$name/" >> windows.fa
  done
  samtools faidx windows.fa
fi
for RUN in $RUNS; do
  if [ ! -s ${RUN}.bam ]; then
    fastq-dump -X 40000000 --stdout $RUN \
      | minimap2 -ax sr windows.fa - 2>/dev/null \
      | samtools view -b -F 4 -q 20 - | samtools sort -o ${RUN}.bam
    samtools index ${RUN}.bam
  fi
done
python3 - << 'PYEOF'
import subprocess, re
sites = [("CRYBB2", 50001, "stop codon"), ("CRYGB", 50001, "stop codon"),
         ("BFSP1", 50001, "frameshift region"), ("BFSP1", 50394, "stop codon"),
         ("BFSP2", 50001, "frameshift + stop")]
for run in ["SRR13005646", "SRR39498552"]:
    print(f"\n=== {run}" + ("  (second individual)" if run == "SRR13005646" else "  (the assembly's own animal)"))
    for win, pos, what in sites:
        lo, hi = pos - 60, pos + 60
        out = subprocess.run(["samtools", "mpileup", "-f", "windows.fa", "-r", f"{win}:{lo}-{hi}", f"{run}.bam"],
                             capture_output=True, text=True).stdout.splitlines()
        depths, indels, mism = [], [], []
        for l in out:
            f = l.split("\t")
            p, d, bases = int(f[1]), int(f[3]), f[4]
            depths.append(d)
            n_ind, n_mis, i = 0, 0, 0
            while i < len(bases):
                c = bases[i]
                if c == "^": i += 2; continue
                if c in "+-":
                    m = re.match(r"[+-](\d+)", bases[i:]); n_ind += 1
                    i += len(m.group(0)) + int(m.group(1)); continue
                if c in "ACGTNacgtn": n_mis += 1
                i += 1
            if n_ind: indels.append(f"{p}({n_ind}/{d})")
            if n_mis and d and n_mis / d >= 0.5: mism.append(f"{p}({n_mis}/{d})")
        cov = f"depth {min(depths) if depths else 0}-{max(depths) if depths else 0}" if depths else "NO READS"
        verdict = "reads AGREE with the assembly" if depths and min(depths) >= 2 and not indels and not mism else "CHECK"
        print(f"  {win:7} {what:18} {cov:14}  indels in reads: {', '.join(indels) or 'none'}   "
              f"majority mismatches: {', '.join(mism) or 'none'}   -> {verdict}")
PYEOF

# ---------------------------------------------------------------- Part C: CRYGC
# CRYGC lies inside the CRYGB window (window positions 59,725-61,399); uses the reads mapped in Part B.
python3 - << 'PYEOF'
import subprocess, re
for run in ["SRR13005646", "SRR39498552"]:
    out = subprocess.run(["samtools", "mpileup", "-f", "windows.fa", "-r", "CRYGB:59700-61450", f"{run}.bam"],
                         capture_output=True, text=True).stdout.splitlines()
    depths, flags = [], []
    for l in out:
        f = l.split("\t"); p, d, b = int(f[1]), int(f[3]), f[4]; depths.append(d)
        ind = len(re.findall(r"[+-]\d+", b))
        if d and ind / d >= 0.3: flags.append(f"indel at {p} ({ind}/{d})")
    low = sum(1 for d in depths if d < 2)
    print(f"CRYGC, {run}: {len(depths)} positions covered, depth {min(depths) if depths else 0}-{max(depths) if depths else 0}, "
          f"{low} positions under 2 reads; {', '.join(flags) or 'no indels in reads'}")
PYEOF
echo "== CRYGC frameshift region, base by base, second individual"
samtools mpileup -f windows.fa -r CRYGB:61330-61385 SRR13005646.bam 2>/dev/null | cut -f2-5
