# river-dolphin-lens-genes

Code for: **Lens genes are lost in the lens-less Ganges river dolphin while light-detection genes are kept**
Muhammad Saeed Iqbal [and co-authors]. [Journal or preprint, year, DOI]

These scripts find gene losses and test for changes in selection in 25 lens, light-detection and hearing
genes, in three river dolphins (*Platanista gangetica*, *Lipotes vexillifer*, *Inia geoffrensis*) and three
marine relatives (*Orcinus orca*, *Tursiops truncatus*, *Phocoena phocoena*). They are the scripts used
for the paper, in the order they were run. No data are stored here: all genomes and reads are public on NCBI.

## Setup

```bash
conda env create -f environment.yml          # main analysis
conda env create -f environment_busco.yml    # BUSCO only (step 2)
conda activate river-dolphin-lens
export LENS_HOME=~/platanista_genomics       # working folder (this is the default)
```

`environment.yml` lists the versions used for the paper. If conda cannot install every pinned version
together on your system, loosen the pins; the versions used are listed below.

| Tool | Version | | Tool | Version |
|---|---|---|---|---|
| BUSCO | 6.1.0 | | SAMtools | 1.19.2 |
| Miniprot | 0.18-r281 | | BCFtools | 1.24 |
| gffread | 0.12.9 | | MACSE | 2.07 |
| minimap2 | 2.31-r1302 | | HyPhy | 2.5.101 |

## Data

| Species | Genome | Reads used to check gene losses |
|---|---|---|
| *Platanista gangetica* | GCA_054084255.1 | SRR13005646 (second individual), SRR39498552 (assembly individual) |
| *Lipotes vexillifer* | GCF_000442215.2 | |
| *Inia geoffrensis* | GCA_036417435.1 | SRR11430609 (second individual), SRR30296906 (assembly project, PacBio HiFi) |
| *Orcinus orca* | GCF_937001465.1 | |
| *Tursiops truncatus* | GCA_011762535.1 | |
| *Phocoena phocoena* | GCF_963924675.1 | |

Genomes are downloaded by step 1. Reads are streamed from the SRA during steps 7 and 8; nothing needs to be downloaded in advance.

## Running the analysis

Run the steps in order from the `scripts/` folder. Each writes into `$LENS_HOME`.

```bash
bash 01_download_genomes.sh
conda activate river-dolphin-busco && bash 02_busco.sh && conda activate river-dolphin-lens
python3 03_lens_pipeline.py fetch
python3 03_lens_pipeline.py map
python3 04_inspect_losses.py
python3 05_check_copies.py
bash 06_rebuild_gamma_crystallins.sh
bash 07_reads_check_platanista.sh
bash 08_reads_check_inia.sh
python3 03_lens_pipeline.py align
python3 03_lens_pipeline.py select
python3 03_lens_pipeline.py report
python3 09_selection_followups.py
python3 10_hearing_checks.py
python3 11_build_supplement.py
python3 12_make_supplement_files.py
```

The slow steps are mapping (step 3 `map`, about an hour), the read checks (steps 7 and 8, a few hours
each, mostly download) and the selection tests (step 3 `select` and step 9, several hours). Steps 7 and 8
can run in the background while the rest continues.

## What each script does

| Step | Script | What it does | In the paper |
|---|---|---|---|
| 1 | `01_download_genomes.sh` | Downloads the six genome assemblies | Table 1 |
| 2 | `02_busco.sh` | Scores genome completeness with BUSCO (cetartiodactyla_odb10) | Table 1 |
| 3 | `03_lens_pipeline.py` | `fetch`: one reference protein per gene from NCBI RefSeq (*O. orca*, then *T. truncatus*, then human; *CDH23* and *TMC1* set explicitly). `map`: maps every reference onto every genome with Miniprot, gives each genome region to one gene only, and calls each gene intact, disrupted (frameshift or premature stop), incomplete (under 80% found) or missing. `align`: codon alignments with MACSE, frameshifts and stops set as gaps. `select`: RELAX for lens and light-detection genes in four designs (all three river dolphins together and each alone); aBSREL on every branch for hearing genes. `report`: Holm correction and summary tables | Figure 3; Tables 4, 5; S2, S3, S4 |
| 4 | `04_inspect_losses.py` | Aligns each disrupted protein to its reference: where the break is, how much protein is lost, and the genome coordinate of the break | Table 3, Figure 4 |
| 5 | `05_check_copies.py` | Lists every match of each broken gene in its genome, to look for an intact second copy, and lists the *P. gangetica* γ-crystallin cluster | Results |
| 6 | `06_rebuild_gamma_crystallins.sh` | Rebuilds *CRYGC* and *CRYGB* in *P. gangetica* from their own loci (the genome-wide model of *CRYGB* ran into the neighbouring *CRYGC*) | Table 3, Table S2 |
| 7 | `07_reads_check_platanista.sh` | Maps reads from a second individual and from the assembly animal around each *P. gangetica* break and checks they agree with the assembly | Table 3, Table S5 |
| 8 | `08_reads_check_inia.sh` | The same for *CRYBA1* in *I. geoffrensis*, with Illumina and PacBio HiFi reads | Table 3, Table S5 |
| 9 | `09_selection_followups.py` | Reruns *CRYBA1*, *BFSP2* and all *CDH23* tests on filtered alignments; RELAX for the rebuilt *CRYGC* and *CRYGB* | Tables 4, 5; S3 |
| 10 | `10_hearing_checks.py` | Checks the *CDH23* sequences and signals, and applies Holm correction across all 26 hearing-gene branch tests | Table 5, Table S4 |
| 11 | `11_build_supplement.py` | Writes Tables S1–S5 (CSV) and Data S1 (alignments) | Supplementary material |
| 12 | `12_make_supplement_files.py` | Makes the supplementary PDF and Excel workbook (Tables S1–S6) | Supplementary material |

## Notes

- **Paths.** All scripts read and write under `$LENS_HOME` (default `~/platanista_genomics`): genomes in
  `data/genomes/`, analysis in `Lens_Project/`, read checks in `Lens_reads/` and `CRYGD_reads/`,
  supplementary files in `Supplementary/`.
- **Read checks.** A break is accepted when at least two reads cross it, no read shows an insertion or
  deletion there (for the deeply covered *I. geoffrensis* break, none in 30% or more of reads), and no
  majority of reads changes a stop codon. The *CRYGD* check uses the first 20 million reads of each run;
  the other *P. gangetica* checks use the first 40 million.
- **Table S5** is compiled by hand from the output of steps 7 and 8.
- **Repeat runs.** HyPhy fits models numerically, so repeated runs can give slightly different values in
  the last decimal places.

## License

MIT (see `LICENSE`).

## Contact

Muhammad Saeed Iqbal, [email]
