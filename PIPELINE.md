# CLsat / DarSat analysis pipeline — code deposit

> **Layout note (2026-10-15).** The stage directories `01_discovery_clustering/` ...
> `07_post_catalogue/` referred to below now live under `satpipe/stages/` with the same names,
> as this repository is being turned into the reusable tool `satpipe` (see `README.md`,
> `HANDOFF.md`). The code exactly as deposited for the manuscript is tag `manuscript-deposit`.
> Superseded code (`superseded_blastn/`, `expand_clsats.sh`) is in `legacy/`.

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
07_post_catalogue/             Singleton screen (L541), candidates C1-C7, IMO-TaqI test, Figure 2
viewer/                        Pointer to the viewer repository; data-product provenance table
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

Step 5 is `Step5.sh` (plus-strand lexicographically-minimal rotation; the reverse complement is
deliberately not considered, matching `vsearch --strand plus` in Steps 6 and 8 — this is what
left L541 unclustered). Its core is:

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
| `build_master.sh` | Union the four families' expanded intervals (`out_CLsat<N>.expanded_merged.bed`, or `EXP_CLSAT<N>=<path>`), merge per species with family provenance (`bedtools merge -c 4,5 -o distinct`) → `MASTER.merged.bed` |

> **Input provenance.** The published MASTER was built from the `cli_v6` re-expansion of CLsat1
> (`v6full_CLsat1.expanded_merged.bed`) and the Jun-12 expansions of CLsat2–4
> (`out_CLsat{2,3,4}_v4.expanded_merged.bed`, whose `MIN_ALN` was 292/290/232 bp against
> cli_v6's 293/291/233 — a 1 bp difference). Until 2026-10-15 `build_master.sh` hard-coded the
> `_v4` name for all four families, which for CLsat1 is the truncated Jun-14 run.
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
consensuses, Δ≥10 bits = confident; bitscores compared as floats since 2026-10-15),
`terminal_monomers.sh` (rewritten 2026-10-15: window vs doubled consensus, strand recorded,
phase read at the array edge — the earlier version reported the best-scoring monomer in the
window, not the terminal one), `boundary_phase.sh`,
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
bedtools v2.31.1 · GNU parallel · Python 3.12.11 (numpy 2.3.2) ·
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

| | De novo | Pipeline output |
|---|---|---|
| Raw MASTER loci | 540 | 540 |
| Segmented loci | 626 | 626 |
| Total span | 68.3995 Mb | 68.4212 Mb |
| CLsat2 / CLsat4 loci | 84 / 15 | 84 / 15 |
| CLsat1 / CLsat3 loci | 429 / 98 | 430 / 97 |

Span agreement: **99.968%**.

> **The manuscript reports 627 segmented / 541 MASTER loci.** The extra locus is **L541**
> (mix, `JAWWNH010000226.1:1–43,510`), which this pipeline does **not** and cannot produce.
> Its TRF monomer is reported in the opposite orientation to the other mix CLsat3 monomers,
> and Step 6 runs `vsearch --strand plus`, so it stayed a singleton locus and never reached
> the cross-assembly pass that generates seeds (Step 8). It was recovered afterwards by
> screening all 979 singleton loci against the four CLsat consensuses — 14 scored above the
> 50-bit floor, seven of which already lay inside catalogued loci — and added by hand.
> Everything above measures the code; the published catalogue is the code's output plus
> that one locus. The seven further singleton loci that the same screen found *outside* the
> catalogue are reported in the manuscript but deliberately excluded from it.

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
   UPGMA tree inspection. It is deposited as an explicit table (above, and as
   `02_seed_expansion_scan/cluster_family.tsv`) so the result is
   regenerable; it is not claimed to be an automatic consequence of a threshold.

### Reproducing this

Steps 3–10 run in minutes. Step 11 (expansion) and the 2nd scan are the expensive
stages — roughly 1.5 h and 3 h respectively at `--jobs 16` on a 64-core host. The
85 Mb rescan requires the bisection cascade described above.

---

## Changes after the reproduction report (2026-10-15 audit)

The reproduction report above describes the code as deposited on 2026-08-04. The following fixes
were applied afterwards; those marked **re-run** can change numbers and have not been re-run
against the genomes.

| Script | Change | Effect on a re-run |
|---|---|---|
| `01/Step4.sh` | `asorti(start, idx)` sorted index *strings* ("1","10","11",…,"2"), so the overlap merge saw starts out of order for motifs with ≥10 TRF hits; now visits intervals by start value. `trf` no longer aborts `set -e`; outputs are truncated before appending. | **re-run**: coverage values change for such motifs, so the ≤50 bp simple-repeat filter may admit or drop a few motifs |
| `01/Step3.sh`, `01/Step7.sh` | `grep -c … \|\| echo 0` printed `0\n0` on no match | none |
| `01/Step6.sh` | cluster files renamed by vsearch's own number, not lexical order | file names only |
| `02/clisat_expand_then_scan_merged.sh` | rows whose scaffold begins with `Chr` are no longer discarded; every locus that cannot be processed is written to `OUT.dropped.tsv` and the seed count is reconciled; header, defaults and `--help` corrected; Phase-2/3 'signal' test documented as weaker than Phase 1 | none on the published assemblies |
| `03/build_master.sh` | reads `out_CLsat<N>.expanded_merged.bed` (overridable), `LC_ALL=C` | see provenance note above |
| `03/arbitration_patched.sh`, `06/rearbitration/arbitration_540.sh` | branch 3 writes Δ = seed − strongest other (same sign as branches 1–2; contradiction when Δ ≤ −20) and visits families in sorted order; D3b uses `(pos−1)//BIN` like D3; short segments are written to `MASTER.loci.segmented.bed.short_segments.bed` (default `SEG_SHORT=drop` keeps the published behaviour; `SEG_SHORT=merge` folds them into the longer neighbour); `LC_ALL=C` | Δ sign only, unless a scan window starts at a multiple of 1000 bp (**re-run** to confirm 626) |
| `04/clean_ortho_pairs.py`, `04/rebuild_ortho_segmented.py` | overlap guard converts the 1-based flank to 0-based before the interval test (1 bp) | **re-run**: a flank starting on the last base of a neighbouring locus is now rejected |
| `04/` | `disc_ortho.py`, `flank_classify.py`, `ortho_matrix.py` moved to `superseded_blastn/consumers/` | none |
| `05/assign_subfamily.sh`, `06/subfamily_rerun/run_assign.sh` | bitscores compared as floats (integer truncation let L299, L327, L396, L407 be called with negative Δ) | **re-run**: those four low-confidence calls flip to the higher-scoring sub-subfamily |
| `05/terminal_monomers.sh`, `build_subfam_term_json.py` | rewritten (see §05); adds `strand` and `edge_gap` | **re-run**: all terminal-monomer values; the viewer's old file is withdrawn |
| `07/imotaq_compare.sh` | DarSat FASTA resolved before `cd`; second argument | none |
| `07/README.md`, `singleton_screen.sh` | singleton-screen arithmetic reconciled (6 inside + L541 + 7 outside = 14) | none |
| deposit | `Step5.sh`, `Step8_cross_cluster.sh`, `Step9_family_fasta.sh`, `cluster_family.tsv`, `consensus/CLsat{1..4}.cons.fa`, `genomes.example.tsv`, `cli_v6.sh`/`cli.sh` links | — |

Still not deposited: `CLsat_cons.rot` (Step 8 references) and the sub-subfamily consensus FASTA
(`subx/`); `sub44_4_4_prefix.fa` in the 85 Mb recovery was made by hand.

## Viewer

The viewer is maintained in [Toki-bio/CLsat-viewer](https://github.com/Toki-bio/CLsat-viewer)
(live at <https://toki-bio.github.io/CLsat-viewer/>); `viewer/README.md` here maps each data
product it reads to the script that produces it. It draws per-window bitscore tracks for all four
subfamilies along each array, so composite loci and subfamily transitions are visible directly
rather than inferred from the arbitration table.
