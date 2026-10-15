# Provenance

This repository was the code deposit for the *Darevskia* CLsat/DarSat manuscript. On
2026-10-15 it was turned, in place, into the reusable tool (`satpipe`, working name). The deposit
as cited by the manuscript is tag **`manuscript-deposit`** (= `main` at `80dd458`, 2026-08-04
scripts + 2026-10-15 documentation). The audit fixes (`345822f`) and the handoff (`4d70b56`) sit
between that tag and the refactor; `git log --follow` on any stage script shows its full history.

## Stage scripts

`satpipe/stages/01_*` ... `07_*` are the former top-level directories, moved with `git mv`
(contents unchanged at the move). The directory names are kept so that line references in the
audit documents and in `PIPELINE.md` still resolve, and because `Step9_family_fasta.sh` and the
stage-06 drivers refer to sibling files by relative path. Three entries are symbolic links:

| link | target |
|---|---|
| `02_seed_expansion_scan/cli_v6.sh` | `clisat_expand_then_scan_merged.sh` |
| `06_catalogue_finalize/cli.sh` | `../02_seed_expansion_scan/clisat_expand_then_scan_merged.sh` |
| `06_catalogue_finalize/cli_v6.sh` | `../02_seed_expansion_scan/clisat_expand_then_scan_merged.sh` |

Moved out of the tool into `legacy/` (HANDOFF.md do-not list: "do not treat ... as current code"),
kept for the record:

- `legacy/superseded_blastn/` (earlier blastn flank-orthology method, 45 pairs instead of 53, plus
  the consumers the audit moved there);
- `legacy/expand_clsats.sh` (serial predecessor of the expansion script, with the `--min-aln-frac` bug);
- `legacy/published_run/` (step 4, 2026-10-15): the one-off drivers of the published catalogue that
  only make sense with the Darevskia files and the paths of that run — the 85 Mb chunk-44 bisection
  cascade (`split_chunk44.sh`, `split_sub44_4.sh`, `bisect_*.sh`, `merge_chunk44_final.py`,
  `proper_verify.py`), the candidate rescans of the recovered locus (`run_new_candidates.sh`,
  `scan_newlocus_other_fams.sh`, `append_newlocus_bits.py`), the pass-3 variants
  (`run_pass3_changed.sh`, `run_pass3_delta.sh`, `combine_bits_v2.py`, `run_arbitration_540.sh`) and
  the CLsat1 re-expansion records (`_v6_full_clsat1.sh`, `run_clsat1_full_rescan.sh`). The methods
  they drove (`run_chunked.sh`, `finish_merge85.sh`, `classify_candidates.py`, `run_pass3.sh`,
  `arbitration_540.sh`, `refine_edges.py`) stay in `satpipe/stages/06_catalogue_finalize/`, parameterised.

## Parameterisation (step 4, 2026-10-15)

Every script under `satpipe/stages/` now takes paths from arguments and parameters from
environment variables whose defaults are the published values. The usage header of each script
documents its arguments. Three tables replace the former hard-coded `declare -A` maps:

| table | columns | used by |
|---|---|---|
| `genomes.tsv` | species code, genome FASTA (indexed) | expansion, rescan, orthology, sub-subfamily, terminal, finalisation, post |
| `families.tsv` | family, MONO, consensus FASTA | build_locus_master, terminal_monomers, boundary_phase, build_subfam_term_json, rebuild_ortho_segmented, clean_ortho_pairs, run_pass3, refine_edges, classify_candidates, candidate_* |
| `members.tsv` | family, sub-subfamily | assign_subfamily, run_assign, candidate_subsubfamily |
| `arbitration_inputs.tsv` | family, MONO, rescan bits TSV, seed BED | arbitration_patched (and its wrapper arbitration_540) |
| `expansions.tsv` | family, expanded_merged BED | build_master, clean_ortho_pairs |

Equivalence with the original scripts was tested on a synthetic project (`tests/toy/`), see
HANDOFF.md progress log.

## Reference data

| file | role |
|---|---|
| `satpipe/stages/02_seed_expansion_scan/consensus/CLsat{1,2,3,4}.cons.fa` | family consensuses (146, 145, 146, 116 bp); `configs/darevskia_clsat.yaml` points here |
| `satpipe/stages/02_seed_expansion_scan/cluster_family.tsv` | curated Step-9 cluster -> family table (Step 9 reads it from its own directory by default) |
| `satpipe/stages/02_seed_expansion_scan/genomes.example.tsv` | format of the `genomes.tsv` the expansion script reads |
| `ref/darevskia/manual_loci.bed` | written 2026-10-15: L541, mix `JAWWNH010000226.1:1-43,510` (1-based inclusive in PIPELINE.md -> BED `0 43510`), family CLsat3 in column 7; the hand-added 627th locus |

Still only on the analysis server (not in any repository): `CLsat_cons.rot` (Step-8 rotated
references), the sub-subfamily consensus directory `subx/`, the six TRF `.dat` files, the genomes.
`configs/darevskia_clsat.yaml` points at their recorded server paths; they were not verified from
the session that wrote it (server unreachable), so run `satpipe check` there first.

## Viewer (step 8, 2026-10-15)

`satpipe/viewer_html/index.html` is the CLsat-viewer page v2.25 (`audit-fixes` @ 109b38c) transformed
by `tools/patch_viewer_config.py` into the config-driven v3.0 page; the same file was committed to
Toki-bio/CLsat-viewer. `satpipe/viewer_lists.py` rebuilds the three inspection lists from the data
files. Verified against the lists that were embedded in v2.25 (`tests/viewer/compare_embedded.py`):
ORTHO_LOCI (91) and FLAGGED_LOCI (35) identical; FULL_TABLE identical in every field except
`both`. **Both-sided:** the manuscript defines high-confidence pairs as "recovered from both flank
sides (18 of 53)"; the shipped `full_flank_ortho_table.tsv` flags 20 loci, but the script that made
it was never deposited and no rule reproduces its flags from the shipped `flank_orthologs.json`
(tested: both own flanks hit; hits in both directions; L+R-sourced support; at 70/80/90 %). The tool
now defines a both-sided pair as one supported by an L-sourced AND an R-sourced hit in either
direction, which on the shipped data gives 6 pairs / 11 loci. The manuscript's "18 of 53" cannot be
reproduced from the deposited files; this is an open point for the owner.

## Documents

- `PIPELINE.md`: the scientific record (reproduction report, changes after it). Unchanged except
  for the layout note at its top.
- `HANDOFF.md`: mission, conventions, target design, work plan, do-not list.
- `docs/config.md`: generated by `satpipe config --schema-doc`.
- The audit this work follows: "CLsat viewer & satellite workflow — code audit"
  (https://claude.ai/code/artifact/91b37065-95a9-49c4-a7b6-8c14d68077db) and "CLsat repos — audit
  fixes applied" (https://claude.ai/code/artifact/976e583b-9cf0-437d-b13a-82059cf3adc6).
