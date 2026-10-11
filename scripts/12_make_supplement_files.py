#!/usr/bin/env python3
"""Step 12. Make the supplementary PDF and Excel workbook (Tables S1-S6) from the CSVs written in step 11.
Output: $LENS_HOME/Supplementary/Supplementary_Material.pdf and Supplementary_Tables.xlsx
Needs reportlab and openpyxl.      Usage: python3 12_make_supplement_files.py
"""
import csv, re, os, glob
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib import colors
from reportlab.lib.units import mm
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.fonts import addMapping
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side

H = os.environ.get("LENS_HOME", os.path.expanduser("~/platanista_genomics"))
U = f"{H}/Supplementary/"
os.chdir(U)
def load(f): return list(csv.DictReader(open(U + f, encoding="utf-8-sig")))
SP = {"Platanista_gangetica": "Platanista gangetica", "Lipotes_vexillifer": "Lipotes vexillifer", "Inia_geoffrensis": "Inia geoffrensis",
      "Orcinus_orca": "Orcinus orca", "Tursiops_truncatus": "Tursiops truncatus", "Phocoena_phocoena": "Phocoena phocoena"}
CAT = {"lens": "lens", "lens-control": "lens control", "photo": "light detection", "photo-control": "light-detection control", "hearing": "hearing"}
def fp(x):
    if x in ("", None): return "–"
    if x.strip().startswith("<"): return x.strip()
    v = float(x); return "< 0.001" if 0 < v < 0.0001 else (f"{v:.3g}" if v < 1 else "1")
def fk(x): return "–" if x in ("", None) else f"{float(x):.2f}"
def it(s): return f"<i>{s}</i>"
SPECIES_PAT = re.compile(r"(Platanista gangetica|Lipotes vexillifer|Inia geoffrensis|Orcinus orca|Tursiops truncatus|Phocoena phocoena|Homo sapiens|P\. gangetica|L\. vexillifer|I\. geoffrensis|O\. orca|T\. truncatus|P\. phocoena)")
def ital_sp(s): return SPECIES_PAT.sub(lambda m: it(m.group(1)), s)

# ---------- clean tables: each column is (header, [cells]); cells hold PDF markup; plain text derived by stripping tags
tables = []
s1 = load("S1_Reference_proteins.csv")
fix = {"TMC1": "transmembrane channel-like protein 1", "CDH23": "cadherin-23",
       "CRYGA": "gamma-crystallin A (not analysed: no separate locus in cetaceans)"}
rows = []
for r in s1:
    prot = fix.get(r["gene"], re.sub(r"\s*\[[^\]]*\]\s*$", "", r["title"]))
    rows.append([it(r["gene"]), r["accession"], it(r["organism"]), r["length"], prot])
tables.append(("S1", "Reference proteins", "Reference protein for each gene: accession, species and length.",
               ["Gene", "Accession", "Species", "Length (aa)", "Protein"], rows, [22, 32, 38, 20, 130]))

s2 = load("S2_Gene_status.csv"); rows = []
for r in s2:
    dash = lambda v: "–" if v in ("-", "") else v
    num = lambda v: "–" if v in ("-", "") else f"{int(v):,}"
    rows.append([it(r["gene"]), CAT[r["category"]], it(SP[r["species"]]), r["status"].lower(), dash(r["scaffold"]), num(r["start"]), num(r["end"]),
                 dash(r["strand"]), ("–" if r["identity"] in ("-", "") else f"{float(r['identity']):.3f}"), ("–" if r["coverage"] in ("-", "") else f"{float(r['coverage']):.2f}"), dash(r["frameshifts"]), dash(r["stop_codons"]), dash(r["exons"])])
tables.append(("S2", "Gene status", "Status of every gene in all six genomes, with coordinates, identity, coverage, frameshifts and stop codons.",
               ["Gene", "Category", "Species", "Status", "Scaffold", "Start", "End", "Strand", "Identity", "Coverage", "Frameshifts", "Stop codons", "Exons"],
               rows, [16, 22, 31, 18, 30, 21, 21, 14, 15, 17, 21, 17, 12]))

s3 = load("S3_RELAX.csv"); rows = []
for r in s3:
    tb = r["test branch"]
    rows.append([it(r["gene"]), CAT[r["category"]], tb if tb == "all three" else it(tb), fk(r["K"]), fp(r["p"]), fp(r["Holm p"]),
                 fk(r["K after filtering"]), fp(r["p after filtering"])])
tables.append(("S3", "All RELAX tests", "All RELAX tests (21 genes, four designs): K, p, Holm-corrected p and filtered reruns. K < 1, relaxed; K > 1, intensified.",
               ["Gene", "Category", "River dolphin branch tested", "K", "p", "Holm p", "K after filtering", "p after filtering"],
               rows, [22, 30, 38, 18, 18, 18, 26, 26]))

s4 = load("S4_aBSREL.csv"); rows = []
hdr4 = list(s4[0].keys())
for r in s4:
    b = r["branch"]
    m = re.match(r"ancestor of (.+)", b)
    if m:
        parts = [p.strip() for p in m.group(1).split("+")]
        b = "Ancestor of " + (" and ".join([", ".join(parts[:-1]), parts[-1]]) if len(parts) > 1 else parts[0])
    rows.append([it(r["gene"]), ital_sp(b)] + [fp(r[h]) if r[h] != "" else "–" for h in hdr4[2:]])
tables.append(("S4", "All aBSREL tests", "All aBSREL branch tests for the three hearing genes, with per-gene and 26-test corrections, before and after filtering.",
               ["Gene", "Branch", "Uncorrected p", "Corrected p (per gene)", "Corrected p (all 26 tests)", "After filtering (per gene)", "After filtering (all 26 tests)"],
               rows, [20, 62, 24, 28, 28, 28, 30]))

s5 = load("S5_Read_evidence.csv"); rows = []
for r in s5:
    res = r["result"].replace("no indel in >=30% of reads", "fewer than 30% of reads differ")
    rows.append([it(r["gene"]), it(r["species"]), r["change"], r["genome position of break"], r["run (individual)"], r["reads at break"].replace("-", "–"), res])
tables.append(("S5", "Read evidence", "Raw-read evidence for each confirmed gene loss.",
               ["Gene", "Species", "Change", "Genome position of break", "Read run (individual)", "Reads at break", "Result"],
               rows, [16, 24, 34, 48, 52, 18, 50]))

rows = [[it("Platanista gangetica"), "SRR13005646", "Illumina", "40 million", "second individual", "PRJNA675309", "SAMN16703806", "HRRL_PG_001"],
        [it("Platanista gangetica"), "SRR39498552", "Illumina", "40 million", "assembly individual", "PRJNA1206299", "SAMN46076716", "GRD9"],
        [it("Inia geoffrensis"), "SRR11430609", "Illumina", "60 million", "second individual", "PRJNA399465", "SAMN07678125", "IniGeo_1_DISCOVAR"],
        [it("Inia geoffrensis"), "SRR30296906", "PacBio HiFi", "0.4 million", "assembly project", "PRJNA1149706", "SAMN32797734", "VGP sample for mIniGeo1"]]
tables.append(("S6", "Read sets", "Read sets used to check gene losses.",
               ["Species", "Run", "Platform", "Reads used", "Role", "BioProject", "BioSample", "Sample name"],
               rows, [34, 24, 20, 20, 30, 26, 26, 40]))

# ---------- PDF
# DejaVu Sans if available (covers every symbol used); otherwise the built-in Helvetica
_dv = glob.glob("/usr/share/fonts/**/DejaVuSans.ttf", recursive=True)
if _dv:
    _d = os.path.dirname(_dv[0])
    for name, f in [("DV", "DejaVuSans.ttf"), ("DV-B", "DejaVuSans-Bold.ttf"), ("DV-I", "DejaVuSans-Oblique.ttf"), ("DV-BI", "DejaVuSans-BoldOblique.ttf")]:
        pdfmetrics.registerFont(TTFont(name, os.path.join(_d, f)))
    for b, i, n in [(0, 0, "DV"), (1, 0, "DV-B"), (0, 1, "DV-I"), (1, 1, "DV-BI")]: addMapping("DV", b, i, n)
    FN, FB = "DV", "DV-B"
else:
    FN, FB = "Helvetica", "Helvetica-Bold"
cell = ParagraphStyle("c", fontName=FN, fontSize=7.2, leading=8.8)
head = ParagraphStyle("h", parent=cell, fontName=FB)
title = ParagraphStyle("t", fontName=FB, fontSize=15, leading=19, spaceAfter=6)
sub = ParagraphStyle("s", fontName=FN, fontSize=10, leading=13, spaceAfter=14, textColor=colors.HexColor("#444444"))
tt = ParagraphStyle("tt", fontName=FN, fontSize=9.5, leading=12.5, spaceAfter=7)
item = ParagraphStyle("i", fontName=FN, fontSize=9.5, leading=13, spaceAfter=5, leftIndent=0)
W = landscape(A4)[0] - 30 * mm

def footer(c, d):
    c.saveState(); c.setFont(FN, 7.5); c.setFillColor(colors.HexColor("#666666"))
    c.drawRightString(landscape(A4)[0] - 15 * mm, 9 * mm, f"Supplementary material · page {d.page}"); c.restoreState()

story = [Paragraph("Supplementary material", title),
         Paragraph("Lens genes are lost in the lens-less Ganges river dolphin while light-detection genes are kept", sub)]
for code, short, legend, *_ in tables:
    story.append(Paragraph(f"<b>Table {code}.</b> {legend}", item))
story += [Paragraph("<b>Data S1.</b> Codon alignments used in the selection tests (FASTA; separate file Data_S1_alignments.zip).", item),
          Paragraph("<b>Data S2.</b> Analysis scripts, available at [repository link].", item), PageBreak()]
for k, (code, short, legend, header, rows, widths) in enumerate(tables):
    cw = [w * W / sum(widths) for w in widths]
    data = [[Paragraph(h, head) for h in header]] + [[Paragraph(str(c), cell) for c in r] for r in rows]
    t = Table(data, colWidths=cw, repeatRows=1)
    t.setStyle(TableStyle([("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#E8ECF1")),
                           ("GRID", (0, 0), (-1, -1), 0.3, colors.HexColor("#9AA5B1")),
                           ("VALIGN", (0, 0), (-1, -1), "TOP"),
                           ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#F7F8FA")]),
                           ("TOPPADDING", (0, 0), (-1, -1), 2), ("BOTTOMPADDING", (0, 0), (-1, -1), 2), ("LEFTPADDING", (0, 0), (-1, -1), 3), ("RIGHTPADDING", (0, 0), (-1, -1), 3)]))
    story += [Paragraph(f"<b>Table {code}.</b> {legend}", tt), t]
    if k < len(tables) - 1: story.append(PageBreak())
doc = SimpleDocTemplate("Supplementary_Material.pdf", pagesize=landscape(A4), leftMargin=15 * mm, rightMargin=15 * mm,
                        topMargin=14 * mm, bottomMargin=15 * mm, title="Supplementary material", author="Muhammad Saeed Iqbal")
doc.build(story, onFirstPage=footer, onLaterPages=footer)

# ---------- Excel (same cleaned values, markup removed)
strip = lambda s: re.sub(r"<[^>]+>", "", str(s))
wb = Workbook(); ws = wb.active; ws.title = "README"
F = lambda **k: Font(name="Arial", **k)
ws["A1"] = "Supplementary material"; ws["A1"].font = F(bold=True, size=14)
ws["A2"] = "Lens genes are lost in the lens-less Ganges river dolphin while light-detection genes are kept"; ws["A2"].font = F(italic=True, size=10)
r0 = 4
for code, short, legend, *_ in tables:
    ws.cell(r0, 1, f"Table {code}").font = F(bold=True); ws.cell(r0, 2, strip(legend)).font = F(); r0 += 1
ws.cell(r0, 1, "Data S1").font = F(bold=True); ws.cell(r0, 2, "Codon alignments used in the selection tests (FASTA; separate file Data_S1_alignments.zip).").font = F(); r0 += 1
ws.cell(r0, 1, "Data S2").font = F(bold=True); ws.cell(r0, 2, "Analysis scripts, available at [repository link].").font = F()
ws.column_dimensions["A"].width = 12; ws.column_dimensions["B"].width = 120
thin = Side(style="thin", color="9AA5B1")
def val(s):
    s = strip(s).replace(",", "") if re.fullmatch(r"[\d,]+", strip(s)) else strip(s)
    try: return int(s) if re.fullmatch(r"\d+", s) else (float(s) if re.fullmatch(r"\d*\.\d+(e-?\d+)?", s) else s)
    except ValueError: return s
for code, short, legend, header, rows, widths in tables:
    sh = wb.create_sheet(f"{code} {short}"[:31])
    sh.cell(1, 1, f"Table {code}. {strip(legend)}").font = F(bold=True)
    for j, h in enumerate(header, 1):
        c = sh.cell(3, j, h); c.font = F(bold=True); c.fill = PatternFill("solid", fgColor="E8ECF1")
        c.alignment = Alignment(wrap_text=True, vertical="top"); c.border = Border(bottom=thin)
    for i, r in enumerate(rows, 4):
        for j, v in enumerate(r, 1):
            c = sh.cell(i, j, val(v)); c.font = F(italic=str(v).startswith("<i>"))
    for j, h in enumerate(header, 1):
        L = max([len(h)] + [len(strip(r[j - 1])) for r in rows])
        sh.column_dimensions[sh.cell(3, j).column_letter].width = min(48, max(9, L + 2))
    sh.freeze_panes = "A4"
wb.save("Supplementary_Tables.xlsx")
print("built PDF and Excel; tables:", [(t[0], len(t[4])) for t in tables])
