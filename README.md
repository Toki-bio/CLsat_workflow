# CLsat workflow

Analysis pipeline for **CLsat** and **DarSat** satellite DNA across six *Darevskia*
lizard genome assemblies — tandem-repeat discovery, locus boundary resolution,
subfamily arbitration, and cross-assembly flank orthology.

Companion code for *"Locus-resolved comparative analysis of minisatellite-scale
satellite DNA in bisexual and hybrid parthenogenetic Darevskia lizards."*

Interactive viewer: **[toki-bio.github.io/CLsat-viewer](https://toki-bio.github.io/CLsat-viewer/)**

---

## Reproducibility

The pipeline was re-run from scratch in a clean directory, starting only from the raw
TRF output, and compared against the published catalogue at every stage:

| | Reproduced | Pipeline output |
|---|---|---|
| MASTER loci | 540 | 540 |
| **Segmented loci** | **626** | **626** |
| Total span | 68.3995 Mb | 68.4212 Mb |

> **Why the manuscript reports 627 / 541, not 626 / 540.** One further locus — **L541**
> (mix, `JAWWNH010000226.1:1–43,510`; ~298 monomers; CLsat3 at 127.8 bits, 33.6 bits clear
> of the next family) — is **not a product of this pipeline and is not reproduced by it**.
> Its TRF monomer is reported in the orientation opposite to the other mix CLsat3 monomers
> (98.6 % identical as reverse complement, 39.7 % same-strand), and per-assembly clustering
> runs `vsearch --strand plus`, so it could not join a cluster, remained a singleton locus,
> and never entered the cross-assembly pass that generates seeds. It was recovered by a
> retrospective screen of all 979 singleton loci against the four CLsat consensuses
> (Methods §3.1) and added to the catalogue by hand. The table above is therefore the
> correct measure of what the code does; the published catalogue is that plus one
> hand-added locus.

Motif extraction through cross-species clustering is **byte-identical** (including the
centroids file, same md5); seed BEDs reproduce exactly (750 / 288 / 101 / 44, zero
coordinate differences). Span agreement **99.968%**.

Remaining differences — eight loci that differ at boundaries on megabase arrays, and
one scaffold-terminal telomeric region that TRF cannot scan in practical time — are
documented explicitly in [PIPELINE.md](PIPELINE.md#reproduction-report) rather than
smoothed over.

---

## Pipeline

```
TRF (.dat)
   │
   ├─ 01  motif extraction → simple-repeat filter → canonical rotation
   │      → per-species clustering (vsearch) → MAFFT rotation
   │
   ├─ 02  cross-species clustering → CLsat1–4 family assignment
   │      → seed BEDs → boundary expansion + monomer-resolution scan
   │
   ├─ 03  MASTER locus universe → symmetric 4-family rescan
   │      → segmentation + family arbitration
   │
   ├─ 04  flank extraction → classification → cross-assembly orthology
   │
   ├─ 05  sub-subfamily calls, terminal monomers, boundary phase
   │
   └─ 06  85 Mb scaffold recovery, edge refinement, re-arbitration
```

| Directory | Stage |
|---|---|
| `01_discovery_clustering/` | Steps 3–7 — motifs to per-species clusters |
| `02_seed_expansion_scan/` | Steps 8–11 — seeds to expanded, scanned arrays |
| `03_master_arbitration/` | Steps 12–14 — MASTER catalogue and family calls |
| `04_flank_orthology/` | Cross-assembly positional orthology from unique flanks |
| `05_subfamily_terminal_boundary/` | Sub-subfamilies, array ends, phase |
| `06_catalogue_finalize/` | Recovery, refinement, final arbitration |
| `07_post_catalogue/` | Checks run after the catalogue: singleton screen (recovered L541), the seven uncatalogued candidate arrays, the IMO-TaqI identity test, Figure 2 |
| `viewer/` | Pointer to the viewer repository and a map of which script produces each data product |

Full per-script documentation, exact commands, and parameters: **[PIPELINE.md](PIPELINE.md)**

---

## Requirements

TRF · VSEARCH · MAFFT · FASTA (`ssearch36`) · EMBOSS · BLAST+ · samtools · bedtools ·
GNU parallel · Python 3 (numpy) · `gawk` (≥ 4.0, for `PROCINFO["sorted_in"]` and `asorti`)

Exact versions used are listed in [PIPELINE.md](PIPELINE.md#software-versions).

## Notes

- Use `clisat_expand_then_scan_merged.sh` (the `cli_v6` boundary-expansion script).
  Earlier `cli.sh` variants contain a `--min-aln-frac` bug that truncates CLsat1 walks.
- Steps 3–10 run in minutes; boundary expansion and the symmetric rescan are the
  expensive stages (~1.5 h and ~3 h at `--jobs 16` on a 64-core host).
- The CLsat family assignment (Step 9) derives from curated tree inspection and is
  deposited as an explicit cluster→family table (`02_seed_expansion_scan/cluster_family.tsv`),
  not re-derived from a threshold.
- Steps 5, 8 and 9 are scripts (`Step5.sh`, `Step8_cross_cluster.sh`, `Step9_family_fasta.sh`);
  `cli_v6.sh` / `cli.sh` referenced by the stage-06 drivers are links to
  `clisat_expand_then_scan_merged.sh`. The four CLsat consensuses are in
  `02_seed_expansion_scan/consensus/`; `genomes.example.tsv` shows the `genomes.tsv` format.
- Not deposited, because they exist only on the analysis server: the phase-rotated reference
  file `CLsat_cons.rot` used in Step 8 (the rotated form of the four consensuses above) and
  the sub-subfamily consensus FASTA (`subx/`) that `assign_subfamily.sh` reads.
- Every script carries the absolute `/staging/...` paths of the original run; set the path
  variables at the top of each script before running elsewhere.
- 2026-10-15 audit fixes that change numbers on a re-run and must be checked against the
  published catalogue: `Step4.sh` interval sort (coverage filter), `assign_subfamily.sh`
  float comparison (four low-confidence sub-subfamily calls), the D3b bin index (`pos-1`),
  `terminal_monomers.sh` (rewritten; the previous output did not measure terminal monomers).
  See PIPELINE.md "Changes after the reproduction report".

## Citation

Please cite the manuscript. Genome assemblies analysed are public under the accessions
listed in the paper.
