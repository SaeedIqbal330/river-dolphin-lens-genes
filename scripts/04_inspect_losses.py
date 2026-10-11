#!/usr/bin/env python3
"""Step 4. Inspect every disrupted gene (Table 3, Figure 4).

For each gene called DISRUPTED by lens_pipeline.py, the protein model is aligned to its reference
to show where the sequence breaks, how far along the protein the break falls, and the genome
coordinate of the break (used for the raw-read checks in steps 7 and 8).

Run after lens_pipeline.py map.   Usage: python3 04_inspect_losses.py
"""
import os, re
from Bio import Align
from Bio.Align import substitution_matrices

W = os.path.join(os.environ.get("LENS_HOME", os.path.expanduser("~/platanista_genomics")), "Lens_Project")
os.chdir(W)


def read_fasta(p):
    d, n = {}, None
    for l in open(p):
        l = l.strip()
        if l.startswith(">"): n = l[1:].split()[0]; d[n] = []
        elif n: d[n].append(l)
    return {k: "".join(v) for k, v in d.items()}


refs = read_fasta("refs.faa")
al = Align.PairwiseAligner(); al.mode = "global"
al.substitution_matrix = substitution_matrices.load("BLOSUM62"); al.open_gap_score = -10; al.extend_gap_score = -0.5


def cds_segments(sp, gene):
    segs, strand = [], "+"
    for l in open(f"selected_{sp}.gff"):
        f = l.rstrip("\n").split("\t")
        if len(f) > 8 and f[2] == "CDS" and re.search(rf"Parent={re.escape(gene)}(;|$)", f[8]):
            segs.append((f[0], int(f[3]), int(f[4]))); strand = f[6]
    segs.sort(key=lambda s: s[1], reverse=(strand == "-"))
    return segs, strand


def genome_pos(segs, strand, nt_index):
    """Map a 0-based position in the spliced coding sequence to a genome coordinate."""
    for sc, a, b in segs:
        L = b - a + 1
        if nt_index < L:
            return sc, (a + nt_index if strand == "+" else b - nt_index)
        nt_index -= L
    return None, None


cands = [l.rstrip("\n").split("\t") for l in list(open("gene_status.tsv"))[1:]]
cands = [c for c in cands if c[3] == "DISRUPTED" and c[1] != "photo-control"]
for c in cands:
    gene, sp = c[0], c[2]
    prot = read_fasta(f"prot_{sp}.fa").get(gene, ""); ref = refs[gene]
    segs, strand = cds_segments(sp, gene)
    core = prot.rstrip(".*")
    stop = next((i for i, a in enumerate(core) if a in ".*"), None)
    print(f"\n=== {gene}  {sp}   ({len(core)} aa model vs {len(ref)} aa reference; frameshifts={c[10]} stops={c[11]})")
    aln = al.align(ref, core.replace(".", "X").replace("*", "X"))[0]
    r, q = aln[0], aln[1]
    # identity in 15-column windows along the reference, to see where the sequence goes out of frame
    win, cols, marks = 15, list(zip(aln[0], aln[1])), []
    for i in range(0, len(cols), win):
        w = [(x, y) for x, y in cols[i:i + win] if x != "-" and y != "-"]
        pid = 100 * sum(x == y for x, y in w) / len(w) if w else 0
        marks.append("#" if pid >= 70 else ("+" if pid >= 40 else "."))
    print("   identity along the protein (# >=70%, + 40-70%, . <40%):")
    print("   " + "".join(marks))
    if stop is not None:
        sc, gpos = genome_pos(segs, strand, stop * 3)
        print(f"   premature stop at residue {stop+1} of {len(core)} ({100*(stop+1)/len(ref):.0f}% along the reference)  -> genome {sc}:{gpos}")
        print(f"   context: ...{core[max(0, stop-12):stop]}[STOP]{core[stop+1:stop+13]}...")
    # first window where identity drops and stays low = likely frameshift position
    low = next((i for i in range(len(marks)) if marks[i] == "." and all(m != "#" for m in marks[i:i+3])), None)
    if low is not None:
        ref_res = sum(1 for x in r[:low * win] if x != "-")
        q_res = sum(1 for y in q[:low * win] if y != "-")
        sc, gpos = genome_pos(segs, strand, q_res * 3)
        print(f"   sequence drops off near reference residue {ref_res} (model residue {q_res})  -> genome {sc}:{gpos}")
    if segs:
        lo, hi = min(s[1] for s in segs), max(s[2] for s in segs)
        print(f"   exons: {len(segs)}; span {segs[0][0]}:{lo}-{hi} ({(hi-lo+1)/1000:.1f} kb)")
