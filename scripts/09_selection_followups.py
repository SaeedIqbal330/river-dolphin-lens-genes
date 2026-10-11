#!/usr/bin/env python3
"""Step 9. Selection-test follow-ups (Tables 4 and 5, Tables S3 and S4).

 1. Robustness: results significant or close to significant after correction (CRYBA1, BFSP2), and all
    CDH23 tests, are rerun on filtered alignments. Filtering removes every codon column with a gap,
    windows of poor alignment (mean pairwise identity < 50% over 5 codons) and low-complexity
    stretches (12 codons with 3 or fewer distinct amino acids).
 2. CRYGC: RELAX with the rebuilt P. gangetica copy added (six species).
 3. CRYGB: RELAX with the rebuilt P. gangetica copy in place of the genome-wide model.
Run after lens_pipeline.py report and 06_rebuild_gamma_crystallins.sh.
Output: Lens_Project/followup/        Usage: python3 09_selection_followups.py
"""
import os, sys, itertools, subprocess
from concurrent.futures import ThreadPoolExecutor
from Bio.Seq import Seq

import importlib.util
_spec = importlib.util.spec_from_file_location(
    "lens_pipeline", os.path.join(os.path.dirname(os.path.abspath(__file__)), "03_lens_pipeline.py"))
lp = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(lp)   # reuse the pipeline's settings and helpers

os.chdir(lp.W); os.makedirs("followup", exist_ok=True)


def filter_alignment(seqs, win_id=5, min_id=0.5, win_lc=12, max_distinct=3):
    names = list(seqs); L = min(len(s) for s in seqs.values()) // 3
    cod = {n: [seqs[n][3*i:3*i+3].upper() for i in range(L)] for n in names}
    aa = {n: [str(Seq(c).translate()) if ("-" not in c and "N" not in c) else None for c in cod[n]] for n in names}
    keep = [True] * L
    for i in range(L):
        if any(aa[n][i] is None for n in names): keep[i] = False
    pairs = list(itertools.combinations(names, 2))
    def col_id(i):
        v = [aa[a][i] == aa[b][i] for a, b in pairs if aa[a][i] and aa[b][i]]
        return sum(v) / len(v) if v else None
    ids, half = [col_id(i) for i in range(L)], win_id // 2
    for i in range(L):
        w = [x for x in ids[max(0, i-half):i+half+1] if x is not None]
        if keep[i] and w and sum(w) / len(w) < min_id: keep[i] = False
    for n in names:
        for s in range(0, L - win_lc + 1):
            w = [x for x in aa[n][s:s+win_lc] if x]
            if len(w) >= win_lc - 2 and len(set(w)) <= max_distinct:
                for i in range(s, s + win_lc): keep[i] = False
    return {n: "".join(c for c, k in zip(cod[n], keep) if k) for n in names}, sum(keep), L


def run(cmd, log):
    env = dict(os.environ, OMP_NUM_THREADS="3")
    with open(log, "w") as f: subprocess.run(cmd, stdout=f, stderr=subprocess.DEVNULL, env=env)
    return log


def relax_jobs(gene, aln_path, designs, tag):
    present = list(lp.read_fasta(aln_path)); tree = lp.prune(lp.TREE, set(present)); jobs = []
    for d in designs:
        test = [s for s in lp.DESIGNS[d] if s in present]
        ref = [s for s in present if lp.SPECIES[s][1] == "reference"]
        t = f"followup/{gene}_{d}_{tag}.nwk"; open(t, "w").write(lp.newick(tree, test, ref) + ";\n")
        a = f"followup/{gene}_{d}_{tag}.fasta"; open(a, "w").write(open(aln_path).read())
        jobs.append((["hyphy", "relax", "--alignment", a, "--tree", t, "--test", "Test", "--reference", "Reference"],
                     f"followup/{gene}_{d}_{tag}.log"))
    return jobs


def add_rebuilt(gene, cds_file, out_prefix):
    """Insert the rebuilt P. gangetica coding sequence into the gene's sequence set and realign with MACSE."""
    pg = "".join(l.strip() for l in open(cds_file) if not l.startswith(">")).upper()
    if len(pg) % 3 == 0 and pg[-3:] in lp.STOPS: pg = pg[:-3]
    seqs = lp.read_fasta(f"alignments/{gene}_unaligned.fasta"); seqs["Platanista_gangetica"] = pg
    lp.write_fasta(seqs, f"followup/{out_prefix}_unaligned.fasta")
    subprocess.run(["macse", "-prog", "alignSequences", "-seq", f"followup/{out_prefix}_unaligned.fasta",
                    "-out_NT", f"followup/{out_prefix}_NT.fasta", "-out_AA", f"followup/{out_prefix}_AA.fasta"],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    al = {k: lp.mask_stops(v.replace("!", "-")) for k, v in lp.read_fasta(f"followup/{out_prefix}_NT.fasta").items()}
    lp.write_fasta(al, f"followup/{out_prefix}.fasta")
    print(f"   {gene}: aligned {len(al)} species, {len(next(iter(al.values())))//3} codons")
    return f"followup/{out_prefix}.fasta"


DESIGNS = ["all", "plat", "lip", "inia"]
jobs = []

print("== 1. Filtered alignments")
for gene in ["CRYBA1", "BFSP2", "CDH23"]:
    seqs, kept, L = filter_alignment(lp.read_fasta(f"alignments/{gene}.fasta"))
    lp.write_fasta(seqs, f"followup/{gene}_filtered.fasta")
    print(f"   {gene}: kept {kept}/{L} codons ({100*kept/L:.1f}%)")
jobs += relax_jobs("CRYBA1", "followup/CRYBA1_filtered.fasta", DESIGNS, "filt")
jobs += relax_jobs("BFSP2", "followup/BFSP2_filtered.fasta", ["plat"], "filt")
cdh = lp.read_fasta("followup/CDH23_filtered.fasta")
t = "followup/CDH23_absrel_filt.nwk"; open(t, "w").write(lp.newick(lp.prune(lp.TREE, set(cdh))) + ";\n")
jobs.append((["hyphy", "absrel", "--alignment", "followup/CDH23_filtered.fasta", "--tree", t], "followup/CDH23_absrel_filt.log"))

print("\n== 2-3. Rebuilt gamma-crystallins")
jobs += relax_jobs("CRYGC6", add_rebuilt("CRYGC", "crygc_cds.fa", "CRYGC6"), DESIGNS, "full")
jobs += relax_jobs("CRYGB", add_rebuilt("CRYGB", "crygb_cds.fa", "CRYGB_new"), DESIGNS, "new")

print(f"\n== Running {len(jobs)} HyPhy jobs, 4 at a time", flush=True)
with ThreadPoolExecutor(max_workers=4) as ex:
    for log in ex.map(lambda j: run(*j), jobs): print("   done:", log, flush=True)

print("\n== RELAX results")
for gene, designs, tag, label in [("CRYBA1", DESIGNS, "filt", "CRYBA1 filtered"), ("BFSP2", ["plat"], "filt", "BFSP2 filtered"),
                                  ("CRYGC6", DESIGNS, "full", "CRYGC (6 species)"), ("CRYGB", DESIGNS, "new", "CRYGB (rebuilt)")]:
    for d in designs:
        K, p = lp.relax_result(f"followup/{gene}_{d}_{tag}.log")
        print(f"   {label:20} {d:5} " + (f"K={K:6.2f}  p={p:.4f}" if p is not None else "RUN FAILED"))
print("\nCDH23 on the filtered alignment: see step 10 (10_hearing_checks.py)")
