#!/usr/bin/env python3
"""Step 10. Checks on the hearing-gene results (Table 5, Table S4).

 1. CDH23 sequences: length, gaps and identity to T. truncatus in each species, to rule out a
    faulty sequence behind a strong signal.
 2. aBSREL on the filtered CDH23 alignment (from step 9), with the tree showing which internal
    branch each "Node" label refers to.
 3. For each significant CDH23 branch, the rate classes (share of sites and their omega), and the
    codons present only in P. phocoena.
 4. Holm correction across all 26 branch tests of the three hearing genes, on the original
    alignments and with the filtered CDH23 alignment.
Run after 09_selection_followups.py.   Usage: python3 10_hearing_checks.py
"""
import os, sys, json

import importlib.util
_spec = importlib.util.spec_from_file_location(
    "lens_pipeline", os.path.join(os.path.dirname(os.path.abspath(__file__)), "03_lens_pipeline.py"))
lp = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(lp)   # reuse the pipeline's settings and helpers

os.chdir(lp.W)
NAMES = {"Node3": "ancestor of Lipotes + Inia", "Node6": "ancestor of Phocoena + Orcinus + Tursiops",
         "Node8": "ancestor of Orcinus + Tursiops"}

print("== 1. CDH23 sequences: length, gaps and identity to T. truncatus")
al = lp.read_fasta("alignments/CDH23.fasta"); tt = al["Tursiops_truncatus"]
for sp, s in al.items():
    cols = [(a, b) for a, b in zip(s, tt) if a != "-" and b != "-"]
    pid = 100 * sum(a == b for a, b in cols) / len(cols) if cols else 0
    print(f"   {sp:22} ungapped {len(s.replace('-', '')):6d} nt   gaps {100*s.count('-')/len(s):5.1f}%   "
          f"identity to T. truncatus {pid:5.1f}%")

print("\n== 2. aBSREL on the filtered CDH23 alignment (corrected p per branch)")
d = json.load(open("followup/CDH23_filtered.fasta.ABSREL.json"))
for b, v in d["branch attributes"]["0"].items():
    print(f"   {NAMES.get(b, b):42} corrected p = {v.get('Corrected P-value')}")
print("   tree with node names:", d.get("input", {}).get("trees", {}).get("0", "not stored"))

print("\n== 3. Rate classes on the significant CDH23 branches (original alignment)")
d = json.load(open("hyphy/CDH23_absrel.fasta.ABSREL.json"))
for b, v in d["branch attributes"]["0"].items():
    p = v.get("Corrected P-value")
    if p is not None and p < 0.05:
        print(f"   {NAMES.get(b, b)}: corrected p = {p}")
        for w, f in v.get("Rate Distributions", []):
            print(f"      omega {w:.2f} on {100*f:.1f}% of sites")
pp = al["Phocoena_phocoena"]; others = [s for k, s in al.items() if k != "Phocoena_phocoena"]
only = [i // 3 + 1 for i in range(0, len(pp), 3) if pp[i:i+3] != "---" and all(o[i:i+3] == "---" for o in others)]
print(f"   codons present only in P. phocoena: {len(only)} {only}")


def absrel(path):
    out = {}
    for b, v in json.load(open(path))["branch attributes"]["0"].items():
        if v.get("Uncorrected P-value") is not None:
            out[b] = (float(v["Uncorrected P-value"]), float(v["Corrected P-value"]))
    return out


print("\n== 4. Holm correction across all branch tests of the three hearing genes")
orig = {g: absrel(f"hyphy/{g}_absrel.fasta.ABSREL.json") for g in ["SLC26A5", "TMC1", "CDH23"]}
filt = absrel("followup/CDH23_filtered.fasta.ABSREL.json")
for label, use_filtered in [("original alignments", False), ("CDH23 after filtering", True)]:
    keys = [(g, b) for g in orig for b in orig[g]]
    pu = [(filt[b][0] if use_filtered and g == "CDH23" and b in filt else orig[g][b][0]) for g, b in keys]
    adj = lp.holm(pu)
    print(f"   -- {label}: {len(keys)} tests")
    for (g, b), u, a in sorted(zip(keys, pu, adj), key=lambda x: x[1])[:6]:
        print(f"      {g:8} {NAMES.get(b, b):42} uncorrected {u:.3g}   corrected (all tests) {a:.3g}" + ("   < 0.05" if a < 0.05 else ""))
