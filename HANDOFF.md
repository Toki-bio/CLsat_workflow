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
- [ ] Driver: plain Python stage runner vs Snakemake.
- [ ] Whether the viewer stays a single file or splits into `index.html` + `viewer.js` + `viewer_config.json` (single file is a feature for deployment; a separate config is needed either way).
- [ ] Second test organism for the generic run.
- [ ] What to do with L186 in the manuscript (seed family retained; contradiction flagged).
- [ ] Whether `molly/` becomes its own repository now or after submission.
