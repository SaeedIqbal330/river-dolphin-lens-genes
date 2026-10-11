#!/usr/bin/env python3
"""
Lens and sensory-gene pipeline: three river dolphins vs three marine references.
Run the stages in order (each is resumable; finished work is reused):
    python3 03_lens_pipeline.py fetch     # reference proteins from NCBI
    python3 03_lens_pipeline.py map       # Miniprot onto all six genomes + gene status table
    python3 03_lens_pipeline.py align     # MACSE codon alignments, frameshifts/stops masked
    python3 03_lens_pipeline.py select    # RELAX (lens + vision) and aBSREL (hearing), in parallel
    python3 03_lens_pipeline.py report    # summary tables, Holm correction, candidate losses
Output goes to $LENS_HOME/Lens_Project (LENS_HOME defaults to ~/platanista_genomics).
Genomes are expected anywhere under $LENS_HOME/data/genomes (see 01_download_genomes.sh).
"""
import os, re, sys, glob, json, time, subprocess, itertools, collections, urllib.request, urllib.parse
from concurrent.futures import ThreadPoolExecutor
from Bio.Seq import Seq

H = os.environ.get("LENS_HOME", os.path.expanduser("~/platanista_genomics"))
W = f"{H}/Lens_Project"
THREADS, JOBS, THREADS_PER_JOB = "8", 4, "3"

SPECIES = {                       # name: (genome accession, role)
    "Platanista_gangetica": ("GCA_054084255.1", "test"),
    "Lipotes_vexillifer":   ("GCF_000442215.2", "test"),
    "Inia_geoffrensis":     ("GCA_036417435.1", "test"),
    "Orcinus_orca":         ("GCF_937001465.1", "reference"),
    "Tursiops_truncatus":   ("GCA_011762535.1", "reference"),
    "Phocoena_phocoena":    ("GCF_963924675.1", "reference"),
}
RIVER = [s for s, (_, r) in SPECIES.items() if r == "test"]
TREE = ("Platanista_gangetica", (("Lipotes_vexillifer", "Inia_geoffrensis"),
        ("Phocoena_phocoena", ("Orcinus_orca", "Tursiops_truncatus"))))   # published topology (McGowen et al. 2020a)
DESIGNS = {"all": RIVER, "plat": ["Platanista_gangetica"], "lip": ["Lipotes_vexillifer"], "inia": ["Inia_geoffrensis"]}

GENES = {   # symbol: (category, protein name used if the symbol search finds nothing)
    "CRYAA": ("lens", "alpha-crystallin A chain"),     "CRYAB": ("lens-control", "alpha-crystallin B chain"),
    "CRYBA1": ("lens", "beta-crystallin A3"),          "CRYBA4": ("lens", "beta-crystallin A4"),
    "CRYBB1": ("lens", "beta-crystallin B1"),          "CRYBB2": ("lens", "beta-crystallin B2"),
    "CRYBB3": ("lens", "beta-crystallin B3"),          "CRYGA": ("lens", "gamma-crystallin A"),
    "CRYGB": ("lens", "gamma-crystallin B"),           "CRYGC": ("lens", "gamma-crystallin C"),
    "CRYGD": ("lens", "gamma-crystallin D"),           "CRYGS": ("lens", "gamma-crystallin S"),
    "MIP": ("lens", "lens fiber major intrinsic protein"), "LIM2": ("lens", "lens fiber membrane intrinsic protein"),
    "BFSP1": ("lens", "filensin"),                     "BFSP2": ("lens", "phakinin"),
    "RHO": ("photo", "rhodopsin"),                     "GNAT1": ("photo", "guanine nucleotide-binding protein G(t) subunit alpha-1"),
    "CNGA1": ("photo", "cGMP-gated cation channel alpha-1"), "GNAT2": ("photo", "guanine nucleotide-binding protein G(t) subunit alpha-2"),
    "SAG": ("photo", "S-arrestin"),                    "OPN1LW": ("photo", "long-wave-sensitive opsin 1"),
    "OPN1SW": ("photo-control", "short-wave-sensitive opsin 1"),
    "SLC26A5": ("hearing", "prestin"),                 "TMC1": ("hearing", "transmembrane channel-like protein 1"),
    "CDH23": ("hearing", "cadherin-23"),
}
REF_ORGANISMS = ["Orcinus orca", "Tursiops truncatus", "Homo sapiens"]
# For two hearing genes the automatic search does not return the protein used in the paper
# (it picks an incomplete or different isoform), so these full-length references are set explicitly.
REF_OVERRIDE = {"CDH23": ("XP_033697516.1", "Tursiops truncatus"),
                "TMC1":  ("XP_049566736.1", "Orcinus orca")}

# ------------------------------------------------------------------ helpers
def read_fasta(p):
    d, n = {}, None
    for l in open(p):
        l = l.strip()
        if l.startswith(">"): n = l[1:].split()[0]; d[n] = []
        elif n: d[n].append(l)
    return {k: "".join(v) for k, v in d.items()}

def write_fasta(seqs, p):
    open(p, "w").write("".join(f">{k}\n{v}\n" for k, v in seqs.items()))

def genome_path(acc):
    key = acc.split("_")[1].split(".")[0]
    hits = [p for p in glob.glob(f"{H}/data/genomes/**/*.fna", recursive=True) if key in os.path.basename(p)]
    return hits[0] if hits else None

def prune(node, keep):
    if isinstance(node, str): return node if node in keep else None
    kids = [k for k in (prune(c, keep) for c in node) if k is not None]
    return None if not kids else (kids[0] if len(kids) == 1 else tuple(kids))

def newick(node, test=(), reference=()):
    if isinstance(node, str):
        return node + ("{Test}" if node in test else "{Reference}" if node in reference else "")
    return "(" + ",".join(newick(c, test, reference) for c in node) + ")"

def holm(pvals):
    order = sorted(range(len(pvals)), key=lambda i: pvals[i]); m = len(pvals); adj = [0.0] * m; run = 0.0
    for rank, i in enumerate(order):
        run = max(run, min(1.0, (m - rank) * pvals[i])); adj[i] = run
    return adj

STOPS = {"TAA", "TAG", "TGA"}
def mask_stops(seq):
    """Replace internal in-frame stop codons (pseudogenes) with gaps so HyPhy accepts the alignment."""
    cod = [seq[i:i+3] for i in range(0, len(seq) - len(seq) % 3, 3)]
    return "".join("---" if c.upper() in STOPS else c for c in cod)

# ------------------------------------------------------------------ stage 1: fetch
EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
def eu(path, **params):
    url = EU + path + "?" + urllib.parse.urlencode(params); last = None
    for _ in range(4):
        try:
            time.sleep(0.4); return urllib.request.urlopen(url, timeout=60).read().decode()
        except Exception as e:
            last = e; time.sleep(3)
    raise last

def candidates(term):
    ids = json.loads(eu("esearch.fcgi", db="protein", term=term, retmax=60, retmode="json"))["esearchresult"]["idlist"]
    if not ids: return []
    res = json.loads(eu("esummary.fcgi", db="protein", id=",".join(ids), retmode="json"))["result"]
    return [(res[i].get("accessionversion", res[i]["caption"]), res[i]["title"], int(res[i]["slen"])) for i in ids if i in res]

def usable(title, name=None):
    t = title.lower()
    if "low quality" in t or "partial" in t: return False
    if name:
        n = name.lower()
        if not t.startswith(n): return False
        nxt = t[len(n):len(n) + 1]
        return nxt in ("", " ", ",", "[")          # excludes "gamma-crystallin D-like"
    return True

def fetch():
    out = []
    for g, (cat, name) in GENES.items():
        if g in REF_OVERRIDE:
            acc, org = REF_OVERRIDE[g]
            fa = eu("efetch.fcgi", db="protein", id=acc, rettype="fasta", retmode="text").splitlines()
            title = fa[0][1:].split(" ", 1)[1] if " " in fa[0] else g
            seq = "".join(l.strip() for l in fa[1:])
            out.append((g, acc, org, len(seq), title, seq))
            print(f"{g:8} {acc:16} {org:20} {len(seq):5d} aa  {title[:60]}  (set explicitly)", flush=True)
            continue
        chosen = None
        for org in REF_ORGANISMS:
            c = [x for x in candidates(f"{g}[Gene Name] AND {org}[Organism] AND refseq[filter]") if usable(x[1])]
            if not c:
                c = [x for x in candidates(f"\"{name}\"[Title] AND {org}[Organism] AND refseq[filter]") if usable(x[1], name)]
            if c:
                chosen = max(c, key=lambda x: x[2]) + (org,); break
        if not chosen:
            print(f"{g:8} NOT FOUND - add a reference by hand"); continue
        acc, title, L, org = chosen
        seq = "".join(l.strip() for l in eu("efetch.fcgi", db="protein", id=acc, rettype="fasta", retmode="text").splitlines() if not l.startswith(">"))
        out.append((g, acc, org, L, title, seq))
        print(f"{g:8} {acc:16} {org:20} {L:5d} aa  {title[:60]}", flush=True)
    write_fasta({g: s for g, _, _, _, _, s in out}, "refs.faa")
    with open("refs_table.tsv", "w") as f:
        f.write("gene\taccession\torganism\tlength\ttitle\n")
        for g, a, o, L, t, _ in out: f.write(f"{g}\t{a}\t{o}\t{L}\t{t}\n")

# ------------------------------------------------------------------ stage 2: map
def parse_gff(path):
    mrna, cds = {}, collections.defaultdict(list)
    for l in open(path):
        if l.startswith("#"): continue
        f = l.rstrip("\n").split("\t")
        if len(f) < 9: continue
        if f[2] == "mRNA":
            a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            tgt = a.get("Target", "").split()
            if len(tgt) < 3: continue
            mrna[a["ID"]] = dict(id=a["ID"], gene=tgt[0], qs=int(tgt[1]), qe=int(tgt[2]), sc=f[0], start=int(f[3]),
                                 end=int(f[4]), strand=f[6], ident=float(a.get("Identity", 0)),
                                 fs=int(a.get("Frameshift", 0)), stop=int(a.get("StopCodon", 0)), line=l)
        elif f[2] == "CDS":
            p = re.search(r"Parent=([^;]+)", f[8])
            if p: cds[p.group(1)].append(l)
    return mrna, cds

def assign(mrna, ref_len):
    """One gene per genome locus: best match (identity x coverage) wins, so paralogs are not double-counted."""
    alns = []
    for m in mrna.values():
        m["cov"] = (m["qe"] - m["qs"] + 1) / ref_len.get(m["gene"], 1)
        alns.append(m)
    alns.sort(key=lambda m: m["ident"] * min(1.0, m["cov"]), reverse=True)
    taken, chosen = [], {}
    for m in alns:
        if m["gene"] in chosen: continue
        if any(m["sc"] == s and m["start"] <= e and m["end"] >= b for s, b, e in taken): continue
        chosen[m["gene"]] = m; taken.append((m["sc"], m["start"], m["end"]))
    return chosen

def status_of(m):
    if m is None: return "MISSING"
    if m["cov"] < 0.8: return "INCOMPLETE"
    if m["fs"] > 0 or m["stop"] > 0: return "DISRUPTED"
    return "INTACT"

def map_stage():
    refs = read_fasta("refs.faa"); ref_len = {g: len(s) for g, s in refs.items()}
    rows = []
    for sp, (acc, role) in SPECIES.items():
        gp = genome_path(acc)
        if not gp:
            print(f"{sp}: genome {acc} not found under data/genomes - download it first"); continue
        gff = f"map_{sp}.gff"
        log = f"map_{sp}.log"
        if not (os.path.exists(gff) and os.path.getsize(gff) > 0 and os.path.exists(log) and "Real time" in open(log).read()):
            print(f"Mapping {len(refs)} proteins onto {sp} ...", flush=True)
            with open(gff, "w") as o, open(f"map_{sp}.log", "w") as e:
                subprocess.run(["miniprot", "-t", THREADS, "--gff", gp, "refs.faa"], stdout=o, stderr=e, check=True)
        mrna, cds = parse_gff(gff)
        chosen = assign(mrna, ref_len)
        sel = f"selected_{sp}.gff"
        with open(sel, "w") as o:
            for g, m in chosen.items():
                o.write(m["line"].replace(f"ID={m['id']}", f"ID={g}", 1))
                for c in cds[m["id"]]: o.write(c.replace(f"Parent={m['id']}", f"Parent={g}"))
        subprocess.run(["gffread", sel, "-g", gp, "-x", f"cds_{sp}.fa", "-y", f"prot_{sp}.fa"], check=True)
        for g in GENES:
            m = chosen.get(g)
            rows.append([g, GENES[g][0], sp, status_of(m)] + ([m["sc"], m["start"], m["end"], m["strand"], f"{m['ident']:.3f}",
                        f"{m['cov']:.2f}", m["fs"], m["stop"], len(cds[m["id"]])] if m else ["-"] * 9))
    with open("gene_status.tsv", "w") as f:
        f.write("gene\tcategory\tspecies\tstatus\tscaffold\tstart\tend\tstrand\tidentity\tcoverage\tframeshifts\tstop_codons\texons\n")
        for r in rows: f.write("\t".join(map(str, r)) + "\n")
    short = {"INTACT": "ok", "DISRUPTED": "LOSS?", "INCOMPLETE": "part", "MISSING": "--"}
    st = {(r[0], r[2]): short[r[3]] for r in rows}
    sps = list(SPECIES)
    print("\n" + f"{'gene':8} {'cat':13} " + " ".join(f"{s.split('_')[0][:10]:>10}" for s in sps))
    for g in GENES:
        print(f"{g:8} {GENES[g][0]:13} " + " ".join(f"{st.get((g, s), '?'):>10}" for s in sps))
    print("\nok = intact; LOSS? = frameshift or stop codon (verify); part = <80% of the protein found; -- = not found")

# ------------------------------------------------------------------ stage 3: align
def align_stage():
    status = collections.defaultdict(dict)
    for l in list(open("gene_status.tsv"))[1:]:
        f = l.rstrip("\n").split("\t"); status[f[0]][f[2]] = f[3]
    cds = {sp: read_fasta(f"cds_{sp}.fa") for sp in SPECIES if os.path.exists(f"cds_{sp}.fa")}
    os.makedirs("alignments", exist_ok=True)
    for g, (cat, _) in GENES.items():
        if cat == "photo-control": continue
        seqs = {}
        for sp in SPECIES:
            s = cds.get(sp, {}).get(g)
            if not s or status[g].get(sp) == "MISSING": continue
            s = s.upper()
            if len(s) % 3 == 0 and s[-3:] in STOPS: s = s[:-3]
            seqs[sp] = s
        if len(seqs) < 4:
            print(f"{g}: only {len(seqs)} species - skipped"); continue
        final = f"alignments/{g}.fasta"
        if os.path.exists(final): continue
        write_fasta(seqs, f"alignments/{g}_unaligned.fasta")
        subprocess.run(["macse", "-prog", "alignSequences", "-seq", f"alignments/{g}_unaligned.fasta",
                        "-out_NT", f"alignments/{g}_macse_NT.fasta", "-out_AA", f"alignments/{g}_macse_AA.fasta"],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        al = {k: mask_stops(v.replace("!", "-")) for k, v in read_fasta(f"alignments/{g}_macse_NT.fasta").items()}
        write_fasta(al, final)
        print(f"{g}: aligned {len(al)} species, {len(next(iter(al.values())))//3} codons", flush=True)

# ------------------------------------------------------------------ stage 4: select
def run(cmd, log):
    if os.path.exists(log) and ("Likelihood ratio test" in open(log).read() or os.path.exists(cmd[3] + ".ABSREL.json")):
        return log + " (reused)"
    env = dict(os.environ, OMP_NUM_THREADS=THREADS_PER_JOB)
    with open(log, "w") as f:
        subprocess.run(cmd, stdout=f, stderr=subprocess.DEVNULL, env=env)
    return log

def select_stage():
    os.makedirs("hyphy", exist_ok=True); jobs = []
    for g, (cat, _) in GENES.items():
        aln = f"alignments/{g}.fasta"
        if not os.path.exists(aln): continue
        present = list(read_fasta(aln))
        tree = prune(TREE, set(present))
        if cat == "hearing":
            t = f"hyphy/{g}_absrel.nwk"; open(t, "w").write(newick(tree) + ";\n")
            a = f"hyphy/{g}_absrel.fasta"; open(a, "w").write(open(aln).read())
            jobs.append((["hyphy", "absrel", "--alignment", a, "--tree", t], f"hyphy/{g}_absrel.log"))
        elif cat in ("lens", "lens-control", "photo"):
            for d, test in DESIGNS.items():
                test = [s for s in test if s in present]
                ref = [s for s in present if SPECIES[s][1] == "reference"]
                if not test or len(ref) < 2: continue
                t = f"hyphy/{g}_{d}.nwk"; open(t, "w").write(newick(tree, test, ref) + ";\n")
                a = f"hyphy/{g}_{d}.fasta"; open(a, "w").write(open(aln).read())
                jobs.append((["hyphy", "relax", "--alignment", a, "--tree", t, "--test", "Test", "--reference", "Reference"],
                             f"hyphy/{g}_{d}.log"))
    print(f"{len(jobs)} HyPhy runs, {JOBS} at a time", flush=True)
    with ThreadPoolExecutor(max_workers=JOBS) as ex:
        for i, log in enumerate(ex.map(lambda j: run(*j), jobs), 1):
            print(f"  [{i}/{len(jobs)}] {log}", flush=True)

# ------------------------------------------------------------------ stage 5: report
def relax_result(log):
    t = open(log).read()
    Ks = [float(x) for x in re.findall(r"parameter \(K\)\s*=\s*([0-9.]+)", t)]
    ps = re.findall(r"Likelihood ratio test \*\*p =\s*([0-9.eE+-]+)\*\*", t)
    if not Ks or not ps: return None, None
    return (Ks[0] if "Restoring" in t else Ks[-1]), float(ps[-1])

def report_stage():
    rows = []
    for g, (cat, _) in GENES.items():
        for d in DESIGNS:
            log = f"hyphy/{g}_{d}.log"
            if os.path.exists(log):
                K, p = relax_result(log)
                if p is not None: rows.append([g, cat, d, K, p])
    print("RELAX (K < 1 relaxed, K > 1 intensified; Holm across genes within each design)")
    with open("relax_results.tsv", "w") as f:
        f.write("gene\tcategory\tdesign\tK\tp\tp_holm\n")
        for d in DESIGNS:
            sub = [r for r in rows if r[2] == d]
            if not sub: continue
            for r, a in zip(sub, holm([r[4] for r in sub])):
                f.write(f"{r[0]}\t{r[1]}\t{d}\t{r[3]:.3f}\t{r[4]:.4g}\t{a:.4g}\n")
                flag = "  <- significant after Holm" if a < 0.05 else ("  (nominal p < 0.05)" if r[4] < 0.05 else "")
                print(f"  {d:5} {r[0]:8} {r[1]:13} K={r[3]:6.2f}  p={r[4]:.4f}  Holm={a:.4f}{flag}")
    print("\naBSREL (Holm-corrected p per branch, from HyPhy)")
    with open("absrel_results.tsv", "w") as f:
        f.write("gene\tbranch\tcorrected_p\n")
        for g, (cat, _) in GENES.items():
            js = f"hyphy/{g}_absrel.fasta.ABSREL.json"
            if cat != "hearing" or not os.path.exists(js): continue
            for b, v in json.load(open(js))["branch attributes"]["0"].items():
                p = v.get("Corrected P-value"); f.write(f"{g}\t{b}\t{p}\n")
                print(f"  {g:8} {b:24} corrected p = {p}{'   <- significant' if p is not None and p <= 0.05 else ''}")
    print("\nCandidate gene losses (verify each against raw reads before reporting):")
    with open("candidate_losses.tsv", "w") as f:
        for l in list(open("gene_status.tsv"))[1:]:
            x = l.rstrip("\n").split("\t")
            if x[3] in ("DISRUPTED", "INCOMPLETE"):
                f.write(l); print(f"  {x[0]:8} {x[2]:22} {x[3]:10} {x[4]}:{x[5]}-{x[6]}  frameshifts={x[10]} stops={x[11]} coverage={x[9]}")

# ------------------------------------------------------------------ main
if __name__ == "__main__":
    os.makedirs(W, exist_ok=True); os.chdir(W)
    stage = sys.argv[1] if len(sys.argv) > 1 else ""
    {"fetch": fetch, "map": map_stage, "align": align_stage, "select": select_stage, "report": report_stage}.get(
        stage, lambda: print(__doc__))()
