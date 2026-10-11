#!/usr/bin/env python3
"""Step 5. Look for an intact second copy of each broken gene, and list the P. gangetica
gamma-crystallin cluster (Methods: "the P. gangetica gamma-crystallin cluster was also checked by hand").

Part 1 lists every Miniprot match of each broken gene in its genome. An intact second copy would
appear as a match with high identity, coverage above 0.8 and no frameshift or stop.
Part 2 lists every match in the gamma-crystallin cluster, which showed that the genome-wide model
for CRYGB had borrowed its last exon from the neighbouring CRYGC (fixed in step 6).

Run after lens_pipeline.py map.   Usage: python3 05_check_copies.py
"""
import os

W = os.path.join(os.environ.get("LENS_HOME", os.path.expanduser("~/platanista_genomics")), "Lens_Project")
os.chdir(W)

LOSSES = [("Platanista_gangetica", g) for g in ["CRYGD", "BFSP2", "CRYBB2", "CRYGB", "CRYGC", "BFSP1"]] + \
         [("Inia_geoffrensis", "CRYBA1")]

ref_len, n = {}, None
for l in open("refs.faa"):
    l = l.strip()
    if l.startswith(">"): n = l[1:]; ref_len[n] = 0
    else: ref_len[n] += len(l)


def attrs(field):
    return dict(kv.split("=", 1) for kv in field.strip().split(";") if "=" in kv)


print("== Part 1: every match of each broken gene in its genome")
for sp, g in LOSSES:
    print(f"-- {g} ({sp})")
    for l in open(f"map_{sp}.gff"):
        f = l.split("\t")
        if len(f) > 8 and f[2] == "mRNA" and f"Target={g} " in f[8]:
            a = attrs(f[8]); t = a["Target"].split()
            cov = (int(t[2]) - int(t[1]) + 1) / ref_len[g]
            print(f"   {f[0]}:{f[3]}-{f[4]}  identity={a.get('Identity')}  coverage={cov:.2f}  "
                  f"frameshift={a.get('Frameshift', '0')}  stop={a.get('StopCodon', '0')}")

print("\n== Part 2: P. gangetica gamma-crystallin cluster (JBLRPT010000008.1, 30.15-30.35 Mb)")
rows = []
for l in open("map_Platanista_gangetica.gff"):
    f = l.split("\t")
    if len(f) > 8 and f[2] == "mRNA" and f[0] == "JBLRPT010000008.1" and 30150000 < int(f[3]) < 30350000:
        a = attrs(f[8])
        rows.append((int(f[3]), int(f[4]), a["Target"].split()[0], a.get("Identity"), a.get("Rank"),
                     a.get("Frameshift", "0"), a.get("StopCodon", "0")))
for s, e, g, i, r, fs, st in sorted(rows):
    print(f"   {s:>11,} {e:>11,}  {(e-s)/1000:5.1f} kb  {g:6} identity={i}  rank={r}  frameshift={fs}  stop={st}")
