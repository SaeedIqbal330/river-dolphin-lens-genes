#!/usr/bin/env python3
"""Step 11. Build Tables S1-S5 (CSV) and Data S1 (codon alignments, zip) from the analysis output.

 S1 reference proteins | S2 gene status in all six genomes | S3 all RELAX tests |
 S4 all aBSREL tests   | S5 raw-read evidence (compiled from steps 7 and 8)
P. gangetica CRYGB and CRYGC use the rebuilt models (step 6) and their selection tests (step 9).
Run after step 10.  Output: $LENS_HOME/Supplementary/      Usage: python3 11_build_supplement.py
"""
import os, re, csv, glob, json, zipfile

H = os.environ.get("LENS_HOME", os.path.expanduser("~/platanista_genomics"))
W = f"{H}/Lens_Project"; OUT = f"{H}/Supplementary"; os.makedirs(OUT, exist_ok=True)
DESIGNS = {"all": "all three", "plat": "P. gangetica", "lip": "L. vexillifer", "inia": "I. geoffrensis"}
NAMES = {"Node3": "ancestor of Lipotes + Inia", "Node6": "ancestor of Phocoena + Orcinus + Tursiops",
         "Node8": "ancestor of Orcinus + Tursiops"}


def tsv(path):
    rows = list(csv.reader(open(path), delimiter="\t"))
    return rows[0], rows[1:]


def holm(p):
    o = sorted(range(len(p)), key=lambda i: p[i]); m = len(p); a = [0.0] * m; r = 0.0
    for k, i in enumerate(o):
        r = max(r, min(1.0, (m - k) * p[i])); a[i] = r
    return a


def relax_result(log):
    t = open(log).read()
    Ks = [float(x) for x in re.findall(r"parameter \(K\)\s*=\s*([0-9.]+)", t)]
    ps = re.findall(r"Likelihood ratio test \*\*p =\s*([0-9.eE+-]+)\*\*", t)
    if not Ks or not ps: return None, None
    return (Ks[0] if "Restoring" in t else Ks[-1]), float(ps[-1])


def fmt(x): return "< 0.001" if x < 0.001 else f"{x:.3g}"


tables = {}

# ---- S1: reference proteins (CRYGA listed last, marked as not analysed)
h, rows = tsv(f"{W}/refs_table.tsv")
rows = [r for r in rows if r[0] != "CRYGA"] + \
       [r[:4] + [r[4] + " (not analysed: no distinct locus in cetaceans)"] for r in rows if r[0] == "CRYGA"]
tables["S1_Reference_proteins"] = (h, rows)

# ---- S2: gene status, with the rebuilt P. gangetica CRYGB and CRYGC models
def window_row(gff, win_start, gene, ref_len):
    n_exons = sum(1 for x in open(gff) if len(x.split("\t")) > 8 and x.split("\t")[2] == "CDS")
    for l in open(gff):
        f = l.split("\t")
        if len(f) > 8 and f[2] == "mRNA":
            a = dict(kv.split("=", 1) for kv in f[8].strip().split(";") if "=" in kv)
            t = a["Target"].split(); cov = (int(t[2]) - int(t[1]) + 1) / ref_len
            return [gene, "lens", "Platanista_gangetica", "DISRUPTED", "JBLRPT010000008.1", str(win_start + int(f[3]) - 1),
                    str(win_start + int(f[4]) - 1), f[6], a["Identity"], f"{cov:.2f}", a.get("Frameshift", "0"),
                    a.get("StopCodon", "0"), str(n_exons)]

ref_len = {r[0]: int(r[3]) for r in tsv(f"{W}/refs_table.tsv")[1]}
h2, r2 = tsv(f"{W}/gene_status.tsv")
fix = {("CRYGB", "Platanista_gangetica"): window_row(f"{W}/crygb.gff", 30236000, "CRYGB", ref_len["CRYGB"]),
       ("CRYGC", "Platanista_gangetica"): window_row(f"{W}/crygc.gff", 30248000, "CRYGC", ref_len["CRYGC"])}
tables["S2_Gene_status"] = (h2, [fix.get((r[0], r[2]), r) for r in r2])

# ---- S3: all RELAX tests; CRYGC and CRYGB from step 9; Holm recomputed within each design
h, rows = tsv(f"{W}/relax_results.tsv")
rows = [r for r in rows if r[0] not in ("CRYGC", "CRYGB")]
for d in DESIGNS:
    for gene, log in [("CRYGC", f"CRYGC6_{d}_full"), ("CRYGB", f"CRYGB_{d}_new")]:
        K, p = relax_result(f"{W}/followup/{log}.log")
        if p is not None: rows.append([gene, "lens", d, f"{K:.3f}", f"{p:.4g}", ""])
filt = {}
for gene in ["CRYBA1", "BFSP2"]:
    for d in DESIGNS:
        log = f"{W}/followup/{gene}_{d}_filt.log"
        if os.path.exists(log):
            K, p = relax_result(log)
            if p is not None: filt[(gene, d)] = (K, p)
out = []
for d, label in DESIGNS.items():
    sub = [r for r in rows if r[2] == d]
    for r, a in zip(sub, holm([float(r[4]) for r in sub])):
        f = filt.get((r[0], d))
        out.append([r[0], r[1], label, r[3], r[4], f"{a:.4g}", f"{f[0]:.3f}" if f else "", f"{f[1]:.4g}" if f else ""])
tables["S3_RELAX"] = (["gene", "category", "test branch", "K", "p", "Holm p", "K after filtering", "p after filtering"], out)

# ---- S4: all aBSREL branch tests; per-gene (HyPhy) and all-26-test (Holm) corrections, before and after filtering
def absrel(path):
    out = {}
    for b, v in json.load(open(path))["branch attributes"]["0"].items():
        if v.get("Uncorrected P-value") is not None:
            out[b] = (float(v["Uncorrected P-value"]), float(v["Corrected P-value"]))
    return out

orig = {g: absrel(f"{W}/hyphy/{g}_absrel.fasta.ABSREL.json") for g in ["SLC26A5", "TMC1", "CDH23"]}
cfilt = absrel(f"{W}/followup/CDH23_filtered.fasta.ABSREL.json")
keys = [(g, b) for g in orig for b in orig[g]]
pool = dict(zip(keys, holm([orig[g][b][0] for g, b in keys])))
poolf = dict(zip(keys, holm([(cfilt[b][0] if g == "CDH23" and b in cfilt else orig[g][b][0]) for g, b in keys])))
s4 = []
for g, b in keys:
    pu, pc = orig[g][b]
    fu = cfilt.get(b) if g == "CDH23" else None
    s4.append([g, NAMES.get(b, b.replace("_", " ")), fmt(pu), fmt(pc), fmt(pool[(g, b)]),
               fmt(fu[1]) if fu else "", fmt(poolf[(g, b)]) if fu else ""])
tables["S4_aBSREL"] = (["gene", "branch", "uncorrected p", "corrected p (per gene, HyPhy)",
                        "corrected p (all 26 hearing tests, Holm)", "corrected p after filtering (per gene)",
                        "corrected p after filtering (all 26 tests)"], s4)

# ---- S5: raw-read evidence, compiled from the output of steps 7 and 8
tables["S5_Read_evidence"] = (["gene", "species", "change", "genome position of break", "run (individual)", "reads at break", "result"], [
    ["CRYGD", "P. gangetica", "frameshift", "JBLRPT010000008.1:30,263,771", "SRR13005646 (second individual)", "2-3", "agrees with assembly"],
    ["BFSP2", "P. gangetica", "frameshift, then stop", "JBLRPT010000001.1:15,310,985", "SRR13005646 (second individual)", "2-9", "agrees with assembly"],
    ["CRYBB2", "P. gangetica", "stop codon", "JBLRPT010000015.1:73,245,827", "SRR13005646 (second individual)", "4-9", "agrees with assembly"],
    ["CRYGB", "P. gangetica", "stop codon", "JBLRPT010000008.1:30,240,533", "SRR13005646 (second individual)", "2-6", "agrees with assembly"],
    ["CRYGC", "P. gangetica", "frameshift, then stop", "JBLRPT010000008.1:~30,251,890", "SRR13005646 (second individual)", "5-9", "agrees with assembly"],
    ["BFSP1", "P. gangetica", "frameshift region, then stop", "JBLRPT010000016.1:51,264,043 and 51,264,436", "SRR13005646 (second individual)", "2-6", "agrees with assembly"],
    ["CRYBA1", "I. geoffrensis", "frameshift", "CM070917.1:21,129,508", "SRR11430609 (second individual)", "19-25", "agrees with assembly; no indel in >=30% of reads"],
    ["CRYBA1", "I. geoffrensis", "frameshift", "CM070917.1:21,129,508", "SRR30296906 (assembly project, PacBio HiFi)", "3-4", "agrees with assembly"],
])

for name, (head, rows) in tables.items():
    with open(f"{OUT}/{name}.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(head); w.writerows(rows)
    print(f"written: {name}.csv ({len(rows)} rows)")

# ---- Data S1: the final codon alignments used in the selection tests
with zipfile.ZipFile(f"{OUT}/Data_S1_alignments.zip", "w", zipfile.ZIP_DEFLATED) as z:
    for p in sorted(glob.glob(f"{W}/alignments/*.fasta")):
        base = os.path.basename(p)
        if "_unaligned" in base or "_macse_" in base or base in ("CRYGC.fasta", "CRYGB.fasta"): continue
        z.write(p, f"alignments/{base}")
    for p, name in [(f"{W}/followup/CRYGC6.fasta", "CRYGC.fasta"), (f"{W}/followup/CRYGB_new.fasta", "CRYGB.fasta")]:
        z.write(p, f"alignments/{name}")
    for p in sorted(glob.glob(f"{W}/followup/*_filtered.fasta")):
        z.write(p, f"alignments/{os.path.basename(p)}")
    print(f"written: Data_S1_alignments.zip ({len(z.namelist())} alignments)")
