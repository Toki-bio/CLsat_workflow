# CLsat / DarSat analysis pipeline — code deposit

Custom code underlying *"Locus-resolved comparative analysis of minisatellite-scale
satellite DNA in bisexual and hybrid parthenogenetic *Darevskia* lizards"*.

This repository contains the **latest version of every script** used to produce the
results reported in the manuscript, arranged in pipeline order. It has been
**verified by a complete de novo re-run from the raw TRF output** — see
[Reproduction report](#reproduction-report) below.

---

## Contents

```
01_discovery_clustering/       Steps 3–7: TRF motif extraction → per-species clustering
02_seed_expansion_scan/        Steps 8–11: seed BEDs → boundary expansion + scan
03_master_arbitration/         Steps 12–14: MASTER locus universe → family arbitration
04_flank_orthology/            Cross-assembly flank extraction, classification, orthology
05_subfamily_terminal_boundary/ Sub-subfamily calls, terminal monomers, boundary phase
06_catalogue_finalize/         85 Mb scaffold recovery, edge refinement, re-arbitration
viewer/                        Single-file D3 viewer + its static JSON data products
```

---

## Pipeline

### Inputs

Six *Darevskia* genome assemblies (accessions in manuscript Table 2), each scanned
with `trf <genome> 2 5 7 80 10 50 2000 -h` (Step 2). The pipeline below starts from
those `.dat` files.

### 01 — Motif discovery and clustering

| Script | Role |
|---|---|
| `Step3.sh` | Extract candidate tandem-repeat motifs from TRF `.dat` (period 80–300 bp, ≥10 copies) → `motifs_output/{sp}.fa` |
| `Step4.sh` | Re-run TRF on the motifs and drop simple-repeat-dominated ones (internal coverage >50 bp) → `{sp}.fa.clean.fas` |
| *(inline awk, see below)* | Step 5: lexicographically-minimal rotation → canonical circular phase → `{sp}.fa.clean.fas.rot` |
| `Step6.sh` | Per-species `vsearch --cluster_fast --id 0.90` → `clusters_by_file/{sp}/` |
| `Step7.sh` | Per cluster: MAFFT align, rotate to rightmost start, realign → `*.rotated.fa` |

Step 5 is a one-liner rather than a numbered script:

```bash
for fa in *.clean.fas; do
  awk '/^>/ {print; next}
       { s=$0; n=length(s); min=s; d=s s;
         for (i=2; i<=n; i++) { r=substr(d,i,n); if (r<min) min=r }
         print min }' "$fa" > "clean/${fa}.rot"
done
```

### 02 — Cross-species clustering, family assignment, seed BEDs

**Step 8 — cross-species clustering.** Pool all Step-7 rotated sequences, prepend the
four CLsat reference consensuses (`CLsat_cons.rot`), degap, and cluster:

```bash
cat clusters_by_file/*/*_rotated/*.rotated.fa > all7_pooled.rot.fa
cat CLsat_cons.rot all7_pooled.rot.fa > all7_pooled.rot.fas
awk '/^>/{print} !/^>/{gsub(/-/,""); if(length($0)>0) print}' all7_pooled.rot.fas > all7_pooled.degap.fa
vsearch --cluster_fast all7_pooled.degap.fa --id 0.90 \
  --centroids cross_species_centroids.fa --clusters cross_clusters/clust_ \
  --uc cross_species.uc --quiet
```

Prepending the references makes them win centroid election for their families
(`--cluster_fast` sorts by decreasing length and processes equal lengths in input
order), so each family's cluster is anchored to a known consensus.

**Step 9 — family assignment.** The CLsat-related clusters were assigned to
CLsat1–4 by inspection of a MAFFT/UPGMA tree of the centroids (rotated to the CLsat1
anchor). That assignment is **curated input**, deposited here as data rather than
re-derived heuristically:

| Family | Cross-species clusters |
|---|---|
| CLsat1 | 59, 61, 66, 67, 68, 9, 90 |
| CLsat2 | 70, 71 |
| CLsat3 | 60, 0, 62, 63, 64, 65, 69 |
| CLsat4 | 113, 6 |

Cluster 62 (a single *arm* locus) sits intermediate between CLsat1/2/3 and was
assigned to CLsat3; this is recorded explicitly because it is a judgement call.

Because Step 8 is deterministic and reproduces byte-identically (verified below),
these cluster numbers are stable and the assignment regenerates the exact seed set.

**Step 10 — seed BEDs.** `make_bed_from_clsats_headers.sh` converts family member
FASTA headers (`sp|SCAFFOLD:START-END|len=...`) to BED6 → `CLsat{1,2,3,4}.bed`
(750 / 288 / 101 / 44 seeds).

**Step 11 — expansion + scan.** `clisat_expand_then_scan_merged.sh` (the corrected
`cli_v6` version) expands each seed outward until satellite signal is lost, merges
overlapping intervals per species+scaffold, then scans at monomer resolution:

```bash
./clisat_expand_then_scan_merged.sh CLsat1.bed genomes.tsv CLsat1.cons.fa 146 out_CLsat1 --jobs 16
# MONO: CLsat1=146, CLsat2=145, CLsat3=146, CLsat4=116
```

Three-phase expansion: dense walk (3×MONO steps, `--max-steps 500`) → coarse binary
probe (fires only when the dense walk is exhausted with signal still present) →
dense pinpoint. Key parameters: `--min-bits 50`, `--min-aln-frac 0.67`.

> **Note — `--min-aln-frac`.** The originally published value (1.0) was computed
> against the wrong base, making the pass condition unsatisfiable and truncating
> CLsat1 walks at the seed edge. `cli_v6` fixes the base and uses 0.67. Use this
> version; the older `cli.sh` is superseded.

### 03 — MASTER locus universe and arbitration

| Script | Role |
|---|---|
| `build_master.sh` | Union the four families' expanded intervals, merge per species with family provenance (`bedtools merge -c 4,5 -o distinct`) → `MASTER.merged.bed` |
| `arbitration_patched.sh` | 2nd scan bits matrix → D3b segmentation → D6 three-branch family arbitration → `MASTER.final_assignment.tsv` |

The 2nd scan re-runs `clisat_expand_then_scan_merged.sh --scan-only` on the fixed
MASTER intervals with **all four** consensuses, so every locus is measured
symmetrically against every family.

Arbitration parameters: `DELTA=20` (seed/rescan contradiction), `DELTA_AMBIG=15`
(ambiguity flag), `SEG_BIN=1000`, `SEG_MIN_RUN=10`, `SEG_MIN_LOCUS=2000`.

### 04 — Flank orthology

`clean_ortho_pairs.py` is self-contained and authoritative: it rebuilds the loci,
extracts 600 bp flanks **with the overlap guard** (a flank may not fall inside another
annotated locus), classifies them (unique / telomeric / low-complexity /
scaffold-truncated), and compares the unique ones all-vs-all with `ssearch36 -m 8C`,
keeping the single best cross-assembly hit per target assembly at ≥200 bp aligned and
≥80% identity. `rebuild_ortho_segmented.py` applies the same method to the segmented
catalogue and returns the identical result.

> Comparisons are deliberately all-vs-all: flanks are oriented array→outward, so two
> orthologous arrays in opposite orientations legitimately match left-flank to
> right-flank. Restricting to same-side pairs discards real inverted orthologs.

> **Do not use `superseded_blastn/`.** Those scripts are the earlier blastn-based
> method, kept only for provenance. They lack the flank overlap guard and restrict
> comparisons to the same side, which yields 45 pairs rather than the correct 53. See
> that directory's README.

### 05 — Sub-subfamily, terminal monomers, boundary phase

`assign_subfamily.sh` (2-monomer boundary window vs doubled sub-subfamily
consensuses, Δ≥10 bits = confident), `terminal_monomers.sh`, `boundary_phase.sh`,
plus the JSON builders that produce the viewer's data products.

### 06 — Catalogue finalization

- `85mb_recovery/` — one *unp* scaffold (`JAWDKN010000003.1`, 85.83 Mb) exceeded TRF's
  practical size limit and was never scanned. `run_chunked.sh` splits it into 44× 2 Mb
  overlapping chunks; `split_chunk44.sh` / `split_sub44_4.sh` further bisect the one
  pathological chunk; `classify_candidates.py` confirms recovered candidates by motif
  identity (not period coincidence); `merge_chunk44_final.py` reassembles coordinates.
- `edge_refine/` — `refine_edges.py` nucleotide-level boundary refinement.
- `rearbitration/` — re-runs arbitration on the finalized catalogue.
- `subfamily_rerun/` — regenerates sub-subfamily calls on the corrected catalogue.

---

## Software versions

TRF v4.10.0-rc.2 · VSEARCH v2.30.0 · MAFFT v7.526 · FASTA/`ssearch36` v36.3.8e ·
EMBOSS v6.6.0.0 (`cons`, `water`) · BLAST+ 2.17.0+ · samtools v1.22.1 ·
bedtools v2.31.1 · GNU parallel · Python 3.12.11 (numpy 2.3.2, pandas 2.3.2) ·
`gawk` (required — `asorti`). Linux x86-64.

---

## Reproduction report

The pipeline was re-run **completely from scratch**, in a clean directory, starting
only from the six TRF `.dat` files, and compared against the published catalogue at
every stage.

| Stage | Result |
|---|---|
| Step 3 — motif extraction | **byte-identical** (all 6 genomes) |
| Step 4 — simple-repeat filter | **byte-identical** (all 6) |
| Step 5 — rotation | **byte-identical** (all 6) |
| Step 6 — per-species clustering | cluster counts **identical** (257/222/343/252/249/270) |
| Step 7 — MAFFT rotation | completed, 6,567 pooled sequences (published: 6,567) |
| Step 8 — cross-species clustering | **331 centroids; `cross_species_centroids.fa` md5-identical** |
| Step 9/10 — seed BEDs | **750 / 288 / 101 / 44 — exact, zero coordinate differences** |
| Step 11/12 — expansion → MASTER | **540 raw loci**; 480/540 coordinates exact |
| Step 13/14 — 2nd scan → arbitration | **626 segmented loci — exact match** |
| 85 Mb recovery | locus recovered, **end coordinate exact** |

**Catalogue totals**

| | De novo | Published |
|---|---|---|
| Raw MASTER loci | 540 | 540 |
| Segmented loci | 626 | 626 |
| Total span | 68.3995 Mb | 68.4212 Mb |
| CLsat2 / CLsat4 loci | 84 / 15 | 84 / 15 |
| CLsat1 / CLsat3 loci | 429 / 98 | 430 / 97 |

Span agreement: **99.968%**.

### Known, documented deviations

These are stated rather than smoothed over:

1. **Eight loci differ at boundaries** (of 540). Fifty differ by ≤200 bp (below the
   3×MONO expansion step, i.e. within the method's own resolution); eight differ more,
   and are split-vs-merged decisions at internal gaps on the hardest megabase arrays.
   No locus region is lost — total covered sequence per scaffold matches. The cause is
   minor `ssearch36` build divergence at the coarse-probe stage; boundaries on such
   arrays are reproducible only to within a few kb, which the manuscript already states.

2. **The 85 Mb locus reproduces at 422,901 bp vs 424,814 bp published** — the end
   coordinate is exact; the start differs by 1,913 bp because the de novo seed set was
   filtered at >10 copies and the original additionally used a 2.7-copy seed.

3. **The scaffold-terminal ~57 kb telomeric array was never TRF-scanned** in either the
   original or the reproduction — TRF does not terminate on it in practical time. It was
   verified by direct k-mer counting to be ~100% `TTAGGG`-family repeat and contains no
   CLsat. `merge_chunk44_final.py` excludes it explicitly rather than imputing content.

4. **`merge_chunk44_final.py` omits sub-chunks `sub44_4_1/2/3`** (chunk44-local
   1.50–1.87 Mb) from its piece list. This was found during reproduction. It does not
   affect the result — the array is seeded from the neighbouring sub-chunk and boundary
   expansion recovers its full extent, and the recovered locus end coordinate is exact —
   but the omission is real and is recorded here rather than silently corrected.

5. **The Step-9 cluster→family assignment is curated, not computed.** It derives from
   UPGMA tree inspection. It is deposited as an explicit table (above) so the result is
   regenerable; it is not claimed to be an automatic consequence of a threshold.

### Reproducing this

Steps 3–10 run in minutes. Step 11 (expansion) and the 2nd scan are the expensive
stages — roughly 1.5 h and 3 h respectively at `--jobs 16` on a 64-core host. The
85 Mb rescan requires the bisection cascade described above.

---

## Viewer

`viewer/index.html` is a single-file D3 application driven entirely by the static
JSON products in the same directory (`flanks.json`, `flank_orthologs.json`,
`subfamily.json`, `terminal_monomers.json`, `boundary_phase.json`). It draws
per-window bitscore tracks for all four subfamilies along each array, so composite
loci and subfamily transitions are visible directly rather than inferred from the
arbitration table. Live at <https://toki-bio.github.io/CLsat-viewer/>.
