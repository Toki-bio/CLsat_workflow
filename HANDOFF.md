# CLsat workflow → reusable satellite-DNA tool: handoff for Claude Code

2026-10-15 · Sergei Kosushkin (written by the Claude audit session; Claude Code reads this through CLAUDE.md)

## Mission and definition of done

Turn the CLsat/DarSat analysis code (Toki-bio/CLsat\_workflow + Toki-bio/CLsat-viewer) into a tool another lab can run on their own genomes and their own satellite families without reading the Darevskia manuscript: one config file, one driver, stage-level inputs/outputs, and the same viewer served from the run's output directory. The published Darevskia result must remain reproducible from the same tool (that is the regression test).

Done means all of the following hold:

- `satpipe run config.yaml` (name to be chosen) executes stages 01–07 in order on a fresh machine with the listed dependencies, resumes after interruption, and writes a run manifest (tool versions, parameters, input checksums).
- The Darevskia config, pointed at the six assemblies, reproduces the published catalogue to the tolerances already documented in PIPELINE.md (byte-identical through Step 10; 626 segmented loci; span agreement ≥ 99.9 %), with the two curated inputs supplied as files.
- A second, unrelated test case (any genome + any TRF-derived monomer family, e.g. a small public vertebrate assembly) runs end to end and the viewer opens on its output.
- No absolute path, species code, consensus length or threshold is hard-coded in any script.
- A reader can go from README to a first run in under an hour.

## Current state (2026-10-15)

Start from the `audit-fixes` branches, not `main`; they contain the corrected scripts and the honest documentation. Read these three documents first: [code audit](https://claude.ai/code/artifact/91b37065-95a9-49c4-a7b6-8c14d68077db), [fixes applied](https://claude.ai/code/artifact/976e583b-9cf0-437d-b13a-82059cf3adc6), and `CLsat_workflow/PIPELINE.md` (its "Changes after the reproduction report" section lists what changes numbers on a re-run).

| Repository | Branch / commit | State |
| --- | --- | --- |
| [Toki-bio/CLsat\_workflow](https://github.com/Toki-bio/CLsat_workflow) | `audit-fixes` @ `345822f` (main @ `80dd458`) | 56 scripts in `01_`–`07_`, bugs fixed, missing scripts deposited, superseded code quarantined. Every script still carries `/staging/...` paths. |
| [Toki-bio/CLsat-viewer](https://github.com/Toki-bio/CLsat-viewer) | `audit-fixes` @ `20ea392` (main @ `4d8ed45`) | Single-file D3 viewer v2.25, worker-based Smith-Waterman, gzip tracks, 634-locus data (627 + C1–C7). Served from the repo root by GitHub Pages. |

Verified in the audit session: new Smith-Waterman identical to old on 31 test loci; gzip loader; Δ normalisation; Step4 coverage fix on synthetic data; D3b drop/merge modes; terminal-monomer phase arithmetic; syntax of every script. **Not verified** (no genomes or aligners in that environment): any stage run end to end; the viewer in a real browser.

Pending items that need the analysis server (`/staging/...`, conda env `bioinfo`): regenerate `terminal_monomers.json` and `subfamily.json` with the fixed scripts; confirm the seed set 750/288/101/44 after the Step4 fix; deposit `CLsat_cons.rot` and the sub-subfamily consensus FASTA (`subx/`); revoke the PAT embedded in two stale local clones (DEPLOY.md); move `molly/` (Poecilia data) out of the viewer repo.

Owner decisions already taken, do not reopen: `SEG_SHORT=drop` is the published behaviour (merge is opt-in); `MIN_MONOMERS_PER_ARRAY=3` in the viewer vs 4 in searsat16 is intentional and documented; the Step-5 plus-strand canonicalisation is deliberate; branch-3 contradictions keep the seed family.

### Progress log

- 2026-10-15 — **step 3 done** (branch `satpipe`): package `satpipe/` (schema, config loader, `satpipe check`,
  tool registry, stage registry, CLI), `configs/darevskia_clsat.yaml`, `docs/config.md`, `docs/files.md`,
  `tests/` (19 pytest tests, synthetic data). Stage scripts moved to `satpipe/stages/` unchanged; superseded
  code to `legacy/`. Server paths in the Darevskia config are from the project's records and **not yet
  verified** (server unreachable from that session) — run `satpipe check configs/darevskia_clsat.yaml` on the
  server before step 4. `satpipe run` / `viewer` are stubs.

- 2026-10-15 — **step 4 done** (branch `step4`, merged to `main`): every stage script takes its inputs from
  arguments, environment variables (published values as defaults) or three small tables —
  `families.tsv` (family, mono, consensus), `genomes.tsv` (code, FASTA), `members.tsv` (family, sub-subfamily);
  `arbitration_patched.sh` takes `arbitration_inputs.tsv` (family, mono, bits, seed) and handles any number of
  families with window = 2 × MONO. No `/staging`, species code, family name, MONO or window width remains in
  code (only in comments). One-off run records of the published catalogue (85 Mb bisection cascade, pass3
  variants, v6 CLsat1 re-expansion, bits combiners) are in `legacy/published_run/`.
  **Verified** on therioserver with `tests/toy/run_toy_equivalence.sh` (synthetic 2-genome project, stages 03–14,
  orthology, sub-subfamily, terminal, phase): every structural output byte-identical between the original
  scripts (`4d70b56`, paths patched by sed) and the parameterised ones. The only differences are ssearch36 bit
  scores, which the ORIGINAL scripts also change between two runs on identical inputs (shuffled-library
  statistics, up to ~15 % on the toy), and the deltas/flags that threshold them; the harness compares those
  within a noise band (`compare_noise.py`). Not covered by the toy: stage 02 (TRF on genomes), finalisation
  scripts, stage 07 scripts (syntax-checked only).

- 2026-10-15 — owner: "run to the victory, use therioserver instead of DRAGEN". **Steps 5–8 done**, 9–10 in progress:
  - `satpipe run` (satpipe/run.py): stages 00 tables, 02 TRF per scaffold in parallel (given `.dat` files
    symlinked instead), 02b chunked TRF + candidate classification of scaffolds above `large_scaffold_bp`
    → recovered seeds, 03–08, 09 curated-table validation + `cross/cluster_report.tsv`, 10 seeds, 11 expansion,
    12 MASTER, 13 rescan, 14 arbitration, 15 catalogue (+ manual loci, tagged MANUAL, numbered after the
    pipeline loci, scanned for their bits), 16 regions, 17 terminal/phase, 18 sub-subfamily (empty table
    when no consensus set), 19 orthology. DONE markers (input fingerprints + config-section hash), `manifest.json`
    (config, tool versions, code commit, counts, centroid md5). `tests/toy/run_driver_toy.sh` passes against the
    equivalence harness.
  - `satpipe viewer` (satpipe/viewer.py): writes `<run>/viewer/` = data files in the viewer's formats (segmented
    catalogue with ids, gzip tracks, ±5 kb flank background scan, JSON products, scaffold lengths, consensuses),
    `viewer_config.json` (families/monos/colours/consensus/sub-subfamily sequences, species order+labels,
    counts, notes, the three inspection lists from satpipe/viewer_lists.py) and the page. The page
    (`satpipe/viewer_html/index.html` = CLsat-viewer v3.0, made by `tools/patch_viewer_config.py` from v2.25)
    embeds nothing run-specific. Verified in headless Chrome (`tests/viewer/cdp_test.js`) on the shipped
    Darevskia data and on the live site; CLsat-viewer `main` now serves v3.0 with a generated `viewer_config.json`.
    ORTHO_LOCI and FLAGGED_LOCI rebuild identically from the data; the FULL_TABLE "both-sided" flag is now
    *pair supported by an L- and an R-sourced hit* (11 loci vs 20 in the former table whose script was never
    deposited; PROVENANCE.md).
  - Darevskia regression on therioserver (`configs/darevskia_clsat_therioserver.yaml`): six GenBank assemblies
    downloaded, TRF regenerated per scaffold (same TRF version), 85.8 Mb unp scaffold handled by stage 02b.
    Not available there: `subx/` (sub-subfamily consensuses, DRAGEN only) and `CLsat_cons.rot` (derived instead).
    Result: see "Regression result" below once the run has finished.
  - **Step 9 done — molly test** (`tests/molly/run_molly.sh`): three families from the revision session's pooled
    Step-8 families (MolSat1 185 bp / 29 units, MolSat2 238 bp / 4, MolSat3 252 bp / 3), both 2021 haplotype
    assemblies (lat, mex), TRF output reused, curated table = the clusters holding the references. End to end in
    ~25 min at 16 jobs: 87 centroids; seeds 29/4/3; 34 MASTER = 34 segmented loci (28 MolSat1, 3, 3; all
    UNIQUE_SEED_WITH_COMPETITOR); 9 anchorable loci, 0 ortholog pairs; viewer assembled and verified in headless
    Chrome (3 families, 2 assemblies, lists 9/0/34, no script errors). `cross/cluster_report.tsv` shows what the
    curator step is for: cluster 1 (4 mex members, 252 bp, 99.2 % identical to MolSat3) is a second MolSat3
    cluster left unassigned by the automatic table.

## The pipeline as it is

Stages run in this order; each row is the contract a refactor must keep. "Hand" marks a curated input that cannot be derived automatically and must become an explicit, validated file.

| Stage | Script(s) | Input | Output | Parameters now hard-coded | Hand |
| --- | --- | --- | --- | --- | --- |
| 02 TRF scan | (command in PIPELINE.md) | genome FASTA per species | `<sp>.dat` | `trf 2 5 7 80 10 50 2000 -h` | — |
| 03 motif extraction | `01/Step3.sh` | `.dat` files (basename = species code) | `motifs_output/<sp>.fa` (headers \`sp | SCF:START-END | len=.. |
| 04 simple-repeat filter | `01/Step4.sh` | motif FASTA | `.clean.fas` / `.repeats.fas` | TRF coverage ≤ 50 bp | — |
| 05 canonical rotation | `01/Step5.sh` | `.clean.fas` | `clean/<file>.rot` | plus strand only | — |
| 06 per-species clustering | `01/Step6.sh` (run inside `clean/`) | `.rot` | `clusters_by_file/<sp>/<sp>_cluster_N.fa` | vsearch `--id 0.90`, `--strand plus` | — |
| 07 MAFFT rotation | `01/Step7.sh` (run inside each cluster dir) | cluster FASTA | `*_rotated/*.rotated.fa` | rightmost-start rotation | — |
| 08 cross-species clustering | `02/Step8_cross_cluster.sh` | all rotated + `CLsat_cons.rot` | `cross_clusters/clust_N`, centroids | `--id 0.90`; references prepended | `CLsat_cons.rot` not deposited |
| 09 family assignment | `02/Step9_family_fasta.sh` + `cluster_family.tsv` | clusters | `CLsat<N>.members.fa` | — | **yes**: cluster→family table from UPGMA tree inspection |
| 10 seed BEDs | `02/make_bed_from_clsats_headers.sh` | members FASTA | `CLsat<N>.bed` (BED6, 0-based) | gawk 3-arg `match` | — |
| 11 expansion + scan | `02/clisat_expand_then_scan_merged.sh` (= `cli_v6.sh`) | seeds, `genomes.tsv`, consensus, MONO | `out_CLsat<N>.expanded_merged.bed`, `.bits.tsv`, `.dropped.tsv` | MONO per family (146/145/146/116), `--min-bits 50`, `--min-aln-frac 0.67`, 3×MONO step, probe 100×MONO, `--max-steps 500`, `--jobs` | — |
| 12 MASTER | `03/build_master.sh` | four expansion BEDs | `MASTER.merged.bed` (+ provenance) | file names via `EXP_CLSAT<N>` | — |
| 13 symmetric rescan | `cli_v6.sh --scan-only` × 4 consensuses on MASTER | MASTER BED | `MASTER_CLsat<N>.bits.tsv` | window 2×MONO, step MONO | — |
| 14 segmentation + arbitration | `03/arbitration_patched.sh` | MASTER, bits, seeds | `MASTER.loci.segmented.bed`, `MASTER.final_assignment.tsv` | DELTA 20, DELTA\_AMBIG 15, SEG\_BIN 1000, SEG\_MIN\_RUN 10, SEG\_MIN\_LOCUS 2000, window widths 292/290/292/232, `SEG_SHORT` | — |
| flank orthology | `04/rebuild_ortho_segmented.py` (authoritative: `clean_ortho_pairs.py`) | segmented loci, genomes | `flanks.json`, `flank_orthologs.json`, matrix | FLANK 600, GUARD 0, ≥ 200 bp, ≥ 80 % id, ssearch36 | — |
| sub-subfamily, terminal, phase | `05/assign_subfamily.sh`, `terminal_monomers.sh`, `boundary_phase.sh`, `build_subfam_term_json.py` | loci, `subx/` consensuses, `all_regions.fa`/`manifest.tsv` from `06/subfamily_rerun/run_extract.sh` | `subfamily.json`, `terminal_monomers.json` | Δ ≥ 10 confident; 1.5-monomer window, EDGE\_TOL 10 | `subx/` not deposited |
| 06 finalisation | `85mb_recovery/`, `edge_refine/`, `rearbitration/`, `subfamily_rerun/` | — | corrected catalogue | chunk 2 Mb / 50 kb overlap; `EDGE_TOL 50` | `sub44_4_4_prefix.fa` made by hand; L541 added by hand |
| 07 post-catalogue | `singleton_screen.sh`, `candidate_*.py`, `imotaq_compare.sh`, `make_fig2_from_data.py` | catalogue, singletons | C1–C7, figure 2 | 50-bit floor | — |

The viewer (`CLsat-viewer/index.html`) reads: `MASTER.merged.bed`, `MASTER.final_assignment.clean.tsv` (no header, 8 columns), `MASTER.locus_bits.tsv`, `MASTER_CLsat{1..4}.bits.tsv.gz`, `MASTER_flank.bits.tsv.gz`, `flanks.json`, `flank_orthologs.json`, `subfamily.json`, optional `terminal_monomers.json`, `scaffold_lengths.tsv`, `CLsat{1..4}.cons.fa`. It embeds three generated lists (`ORTHO_LOCI`, `FLAGGED_LOCI`, `FULL_TABLE`) and the four consensus strings in the HTML, plus the sub-subfamily cluster sequences (`CLUSTER_SEQS`); those are the parts that are Darevskia-specific.

## Conventions to preserve

- **Coordinates.** TRF `.dat`, Step-3 headers, `samtools faidx` regions, ssearch36 `-m 8C` columns 7–10 and scan window `pos` are 1-based inclusive; every BED and every locus `start/end` in JSON is 0-based half-open. The conversions in `make_bed_from_clsats_headers.sh:25`, `clisat…:263–297,559–579` and `arbitration_patched.sh` D3/D3b (`pos−1`) are the reference. `flanks.json` flank `start/end` are 1-based for L-loci (from `rebuild_ortho_segmented.py`) and 0-based for C-loci (from `candidate_scan.py`): fix the latter when you touch it. The viewer passes BED bounds straight to efetch (1-based), so monomer coordinates it reports are 1-based; keep or fix, but say which.
- **Δ.** Column 6 of the assignment TSV = assigned-family bits − strongest other family, for every flag (since 2026-10-15). Negative means another family scored higher. `UNIQUE_SEED_CONTRADICTION` fires at Δ ≤ −20 and the seed family is still the catalogue call.
- **Locus IDs.** `L<n>` in MASTER order (`"L"NR` over a `LC_ALL=C`-sorted file); segments `L<n>.<k>` with gaps allowed when a short segment was dropped; `C<n>` for uncatalogued candidates. Natural sort in the viewer handles both.
- **Flags.** The arbitration script emits exactly: `UNIQUE_SEED_ONLY`, `UNIQUE_SEED_WITH_COMPETITOR`, `UNIQUE_SEED_CONTRADICTION`, `SEED_CONFLICT_RESOLVED[_AMBIGUOUS]`, `SEGMENT_RESOLVED[_AMBIGUOUS]`, `NO_SEED_RESOLVED`, `NO_SEED_AMBIGUOUS`; the viewer adds `CANDIDATE_UNCATALOGUED`. Keep the viewer's `FLAG_COLORS`/labels/tooltips in step with any addition.
- **Canonicalisation.** Step 5 rotation and vsearch clustering are plus-strand only; this is deliberate and documented (it is what left L541 unclustered). A universal tool should expose `strand: plus|both` and default to the published choice for the Darevskia config.
- **MONO per family** drives window sizes everywhere (3×MONO walk, 2×MONO scan window, D3 widths 2×MONO, viewer `CLSAT_STEP`). It must come from one place.
- **Ortholog pairs** are counted per locus pair (53), not per flank link (76); best hit per target assembly is applied to the matrix, not to the JSON.
- **Published behaviours** that look like bugs but are frozen: short segments dropped (`SEG_SHORT=drop`), coarse-probe phase accepts any ≥ 50-bit hit in a 14.6 kb block (documented in the script header), scaffold-end windows annotated up to 437 bp unscanned. Change only behind options, defaults unchanged.

## Target design

Keep the stage scripts as the units of work (they are tested by the reproduction); wrap them, do not rewrite them. Three additions make the tool reusable:

**1. One config file** (`config.yaml`), the only place paths and parameters live:

```yaml
run_name: darevskia_clsat
out_dir: runs/darevskia_clsat
genomes:                       # species code -> FASTA (code is used in locus names)
  arm: /data/D_armeniaca.fna
  mix: /data/D_mixta.fna
families:                      # one block per satellite family
  CLsat1: {consensus: ref/CLsat1.cons.fa, mono: 146}
  CLsat2: {consensus: ref/CLsat2.cons.fa, mono: 145}
discovery:   {trf: "2 5 7 80 10 50 2000", period: [80, 300], min_copies: 10, simple_repeat_cov: 50, strand: plus, cluster_id: 0.90}
family_assignment: ref/cluster_family.tsv      # curated; validated against the Step-8 clusters
expansion:   {min_bits: 50, min_aln_frac: 0.67, step_mono: 3, win_mono: 3, probe_mono: 100, max_steps: 500, jobs: 16}
arbitration: {delta: 20, delta_ambig: 15, seg_bin: 1000, seg_min_run: 10, seg_min_locus: 2000, short_segments: drop}
orthology:   {flank_bp: 600, guard_bp: 0, min_aln: 200, min_pid: 80}
subfamily:   {consensuses: ref/subx/, confident_delta: 10}
manual_loci: ref/manual_loci.bed               # optional, e.g. L541; recorded in the manifest
```

**2. One driver** (`satpipe`, a Python package with a thin CLI; Snakemake is acceptable if the owner prefers, but a plain stage runner is enough): `satpipe run config.yaml [--from STAGE] [--to STAGE]`, `satpipe check config.yaml` (dependencies, inputs, config schema), `satpipe viewer <run_dir>` (assemble the viewer directory). Each stage is a function that takes the config and the run directory, calls the existing script with explicit arguments and environment variables, writes `<run>/<stage>/DONE` with the input checksums, and is skipped if DONE is newer than its inputs. All stage scripts get their paths from arguments or env, never from constants.

**3. Stage contracts** as files, so partial reuse is possible: a lab that already has seeds starts at stage 11 by placing `seeds/<family>.bed` in the run directory. Document every intermediate file's columns in `docs/files.md` (there are about twenty).

**Viewer generalisation** (`CLsat-viewer`): replace the embedded Darevskia constants with a `viewer_config.json` written by `satpipe viewer`: family names, colours, consensus sequences, MONO per family, species codes → display names, and the three list tables (`ORTHO_LOCI`, `FLAGGED_LOCI`, `FULL_TABLE`) generated from the JSON products instead of pasted into the HTML. `CLUSTER_SEQS`/`CLSAT_TO_CLUSTERS` become optional (`subfamilies` key). The page title and manual text come from the config. The Smith-Waterman, charts and loaders need no change.

**Hand steps become validated inputs**: `cluster_family.tsv` is checked at stage 09 (every listed cluster exists; every family has ≥ 1 cluster; a report lists unassigned clusters with their sizes and best-consensus bits so the curator has something to look at); `manual_loci.bed` entries are tagged `MANUAL` in provenance and listed in the manifest.

## Work plan

Do these in order; each step ends with the check named, and the regression target is always the Darevskia catalogue (626 segmented loci from the pipeline + L541 manual).

1. **Merge `audit-fixes` into `main`** in both repositories (owner approves). Check: Pages serves v2.25; `wc -l MASTER.merged.bed` = 634.
2. **Pending server work** (needs `/staging`): regenerate `subfamily.json` and `terminal_monomers.json`; re-run Steps 4–10 and compare seed BEDs to 750/288/101/44. Check: diff of seed BEDs; if the seed set changed, record it in PIPELINE.md before anything else.
3. **Create the tool repository** (`Toki-bio/satpipe`, or rename `CLsat_workflow`): move stage scripts into `satpipe/stages/`, add `config.yaml` schema, `satpipe check`. Check: `satpipe check` passes on the Darevskia config with all paths present.
4. **Parameterise scripts**: strip every `/staging` constant (112 occurrences — `grep -rn '/staging/'`), every `declare -A` species/consensus table, every literal MONO or window width; feed them from the driver. Keep script names and argument orders. Check: `grep -rn '/staging\|JAWWN\|146\b' satpipe/stages` returns nothing that is not a comment.
5. **Driver, stages 03–14**: wire the existing scripts, DONE markers, manifest. Check: Darevskia run from `.dat` files reproduces `cross_species_centroids.fa` md5 and 540 MASTER loci.
6. **Driver, orthology / subfamily / terminal**: replace `extract_flanks.sh` dependence with `run_extract.sh` logic inside the stage; `subx/` from config. Check: 53 locus pairs, 18 both-sided.
7. **Finalisation as optional stages**: 85 Mb chunked rescan (`edge_refine`, `85mb_recovery`) behind a `large_scaffold_bp` threshold; rearbitration as a re-run of stage 14. Check: L541 and the 85 Mb locus end coordinate reproduce.
8. **Viewer config**: generate `viewer_config.json` and the three lists from the run outputs; load consensus strings and family colours from it. Check: the Darevskia build renders identically to v2.25 (same counts in all three lists; spot-check L1, L186, L541, C1).
9. **Second test case**: pick a small public assembly with a known satellite (document which and why), run end to end with its own config. Check: the viewer opens on its output; README walk-through takes < 1 h.
10. **Packaging**: `environment.yml` (trf, vsearch, mafft, fasta3, emboss, blast, samtools, bedtools, parallel, gawk, python ≥ 3.10), CI that runs `satpipe check` and the unit tests from the audit (Smith-Waterman equivalence, Step4 coverage, D3b modes, terminal-monomer phase) on synthetic data, a `CITATION.cff`, a Zenodo release. Check: a clean container builds the environment and passes CI.

## Environment, tooling and gotchas

- **Where the data is.** Genomes, TRF output, seeds, `subx/`, and all intermediates are on the analysis server under `/staging/Sergei/Lacertidae/Darevskia/Satellites/...` and `/staging/tmp/clsat_v6/`; binaries in the conda env `bioinfo` (`/staging/conda/envs/bioinfo/bin`) and `/staging/miniconda3/bin/trf`. The Claude Code cloud container has none of this: no genomes, no ssearch36/vsearch/mafft/trf. Anything that needs them must run on the server (desktop-app device bridge or the owner running commands) or be mocked with synthetic data.
- **Versions used for the published run**: TRF 4.10.0-rc.2, VSEARCH 2.30.0, MAFFT 7.526, FASTA 36.3.8e, EMBOSS 6.6.0, BLAST+ 2.17.0, samtools 1.22.1, bedtools 2.31.1, Python 3.12.11, numpy 2.3.2. `gawk` is required (`PROCINFO["sorted_in"]`, `asorti`, 3-argument `match`); plain `awk`/`mawk` breaks Steps 4 and 10. The audit container had only mawk; a working gawk 5.2.1 was obtained by unpacking the Ubuntu `.deb` plus `libsigsegv2` — expect to do the same.
- **TRF exit status** is not 0 on success in every build; `Step4.sh` now has `|| true`. The 85 Mb `unp` scaffold cannot be scanned by TRF in one piece (hence `85mb_recovery/`).
- **ssearch36 `-m 8C`** reports one alignment per library sequence per strand; strand is read from subject coordinates (col 9 > col 10 = minus). `searsat16_reference.sh` reads it from the query columns instead — do not copy that.
- **Locale.** Locus numbering depends on `sort` order; every script that sorts now sets `LC_ALL=C`. Keep it in the driver's environment.
- **GitHub access from Claude Code**: both repos are attached with push rights in this session's scope; the git proxy only serves anonymous reads of public repos, and GraphQL / account-wide REST endpoints are blocked (`gh api repos/{owner}/{repo}/...` works, `gh repo list` does not). Shallow clones do not fetch non-default branches: add a fetch refspec for `audit-fixes` or `git fetch origin audit-fixes:refs/remotes/origin/audit-fixes` before relying on tracking.
- **GitHub Pages** does not set `Content-Encoding`, so the viewer inflates `.tsv.gz` with `DecompressionStream`; keep files under 100 MB each; `.nojekyll` must stay.
- **NCBI efetch** is called from the browser with `tool=clsat-viewer&email=...`; sequences for a non-GenBank assembly will need a local FASTA server option (the removed `CONTROL_BASE` stub in `index.html` was that idea) — plan it in the viewer-config step.
- **Owner preferences** (from memory): act rather than ask; verify before claiming; cite line numbers; no narration of mistakes; search the web when versions matter; before any Monsoon (NAU) connection ask him to connect Cisco AnyConnect to vpn.nau.edu — but this project runs on the Tashkent server, not Monsoon.

## Do-not list and open decisions

Do not:

- rewrite stage scripts in another language before the driver reproduces the catalogue with the existing ones — the reproduction report is the only test suite that exists;
- change any default that alters the published numbers (`SEG_SHORT`, probe-phase test, scaffold-end handling, plus-strand canonicalisation, `MIN_MONOMERS_PER_ARRAY`); add options with the published value as default;
- commit the uncompressed bitscore TSVs, `index_*` archives, or `__pycache__`;
- put generated lists back into `index.html` by hand once the viewer reads them from `viewer_config.json`;
- treat `expand_clsats.sh`, anything under `superseded_blastn/`, or `test_msa.py` (deleted) as current code;
- `rm -rf` from command substitutions (blocked by the harness) — list paths explicitly.

Decisions the owner must make before or during the refactor:

- [x] Repository layout — **decided 2026-10-15 by the owner: the tool replaces `CLsat_workflow` in place.** The deposit is frozen as tag `manuscript-deposit` (= former `main` @ `80dd458`); stage scripts moved with `git mv` to `satpipe/stages/`; superseded code in `legacy/`. The GitHub repository name is unchanged (renaming is still open; GitHub redirects a renamed repository).
- [ ] Tool name and licence (none declared today; README says "academic and research use").
- [x] Driver — **decided 2026-10-15 by the owner: plain Python stage runner** (one function per stage calling the existing script; DONE markers with input checksums; manifest). A Snakemake file can be layered on the same stage functions later if a cluster version is ever wanted.
- [ ] Whether the viewer stays a single file or splits into `index.html` + `viewer.js` + `viewer_config.json` (single file is a feature for deployment; a separate config is needed either way).
- [x] Second test organism — **decided 2026-10-15 by the owner: the Amazon molly, *Poecilia formosa*.** The
  CLsat-revision session already ran TRF (`2 5 7 80 10 50 2000 -h`) and workflow Steps 3-8 on the two 2021
  parental-haplotype assemblies, GCA_013036135.2 (*P. latipinna*-derived half, 976 sequences, 674.6 Mb) and
  GCA_013036085.2 (*P. mexicana*-derived half, 1,796 sequences, 745.3 Mb): 406 + 630 motifs, 206 + 355 clean,
  99 + 134 per-half clusters, 45 pooled families with >= 3 units; a recurrent ~186 bp unit in both halves. Data:
  `therioserver:/home/toki/clsat_revision/molly_trf/` (`split/<acc>/*.dat`, `all_trf.tsv`, `wf37/`, `wf37/step8/`);
  family alignments in Toki-bio/CLsat-viewer `molly/alignments/` (commit 4d8ed45). Caveat: those assemblies are
  PacBio CLR. For the satpipe test the molly config needs: the two halves as two "genomes" (codes e.g. `lat`,
  `mex`), the existing `.dat` files as `inputs.trf_dat`, one or more families chosen from the 45 (consensus +
  `mono` from the pooled Step-7 alignment; the ~186 bp family is the obvious first), and a curated
  `family_assignment` table written after inspecting the Step-8 centroids, exactly as for Darevskia.
- [ ] What to do with L186 in the manuscript (seed family retained; contradiction flagged).
- [ ] Whether `molly/` becomes its own repository now or after submission.

## Anchor phasing: DONE 2026-10-16 (request from the CLsat revision session)
Implemented as stage `21_anchor`, `satpipe anchor [--check]` and `satpipe/anchor.py` (commits 5f92a20, 73f6d81); config block `anchor` (family, motif, max_motif_mismatch, min_identity, independent). The nine sub-subfamily consensuses re-rotate to exactly fig1b_rotated.fasta; the Darevskia run passes (4 families in frame). Still to do on the manuscript side: re-deposit `subsubfamily_consensuses.fasta` via `satpipe anchor`, and re-make the molly alignments from anchored consensuses (molly families are unrelated: `anchor.independent: true` keeps each family's own frame).
Original request:
- Anchor phasing must be a mandatory automatic stage: rotate every written consensus/alignment (families, sub-subfamilies, molly) to the configured anchor motif (manuscript Methods 3.4: doubled consensus vs anchor, water, AAAGAAG check) and fail if any output does not start at the anchor.
- Do NOT deposit C:\work\CLsat\manuscript\subsubfamily_consensuses.fasta as is (36 letters out of frame); use fig1b_rotated.fasta or re-rotate by 36.
- Molly alignments in CLsat-viewer molly/alignments/ were made without an anchor step (in phase within families, arbitrary between).

## Regression result (Darevskia, therioserver, 2026-10-16, code f77448a)
Full run from the six GenBank assemblies (TRF regenerated, 85 Mb scaffold via chunked TRF), published cluster table re-derived with tests/darevskia/recurate_from_catalogue.py.
- Seeds CLsat1 771 / CLsat2 288 / CLsat3 101 / CLsat4 44 (published 750/288/101/44; +21 CLsat1 from the 85 Mb scaffold, +2 other).
- MASTER 542 (published 540): +2 = three loci split and one merged at 0.8-7 kb gaps between expansions (unp JAWDKN010000003.1; val JAWWNG010000434.1, 495.1 split; val 030.1 merged). Every locus overlaps a published one in both directions; the only absentee is manual L541, added after MASTER.
- Segmented 627 (published 626 + L541 = 627), catalogue 628 (with the manual locus), ortholog pairs 53 (=), candidates 7 (=), 95 anchorable, 20 flagged, 27 contradiction rows.
- Viewer built by `satpipe viewer`, headless-Chrome test passes (635 loci, 4 families, 3 lists 95/27/628; only favicon 404).
- Not run here: sub-subfamily stage (subx/ only on DRAGEN -> empty table).
