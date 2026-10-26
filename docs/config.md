# satpipe configuration reference

Generated from `satpipe/schema.py` by `satpipe config --schema-doc`; do not edit by hand.

Relative paths are resolved against the directory that holds the config file. Unknown keys are errors. Defaults are the values of the published Darevskia run.

| key | type | required | default | meaning |
|---|---|---|---|---|
| `run_name` | str | yes |  | Name of the run; used in the manifest and in output file names. |
| `out_dir` | path | yes |  | Run directory. Every stage writes to <out_dir>/<stage>/. |
| `genomes` | map -> path | yes |  | Species/assembly code -> genome FASTA. The code is used in every locus name (`code|scaffold:start-end`), so keep it short and free of `|`, `:` and `-`. |
| `inputs.trf_dat` | map -> path |  |  | Code -> TRF `.dat` file for that genome. When given for every genome, stage 02 (TRF scan) is skipped. The published Darevskia run starts here. |
| `inputs.seeds` | map -> path |  |  | Family -> seed BED6 (0-based). When given for every family, stages 03-10 are skipped and the run starts at expansion (stage 11). |
| `families.<key>` | map of blocks | yes | | One block per satellite family. |
| `families.<key>.consensus` | path | yes |  | Single-record FASTA with one monomer of the family consensus. |
| `families.<key>.mono` | int | yes |  | Monomer length in bp. Drives every window size (3xMONO walk, 2xMONO scan window, D3 widths, viewer step). Should equal the consensus length; `satpipe check` warns if it does not. |
| `families.<key>.colour` | str |  |  | Display colour for the viewer (hex). Optional. |
| `references.rotated_consensus` | path |  |  | Phase-rotated family consensuses prepended in Step 8 (`CLsat_cons.rot` in the Darevskia run). Prepending makes each reference win centroid election for its family. If absent the driver derives it by applying the Step-5 rotation to each `families.*.consensus` and records that deviation in the manifest. |
| `discovery.trf` | str |  | `2 5 7 80 10 50 2000` | TRF parameters: match mismatch delta PM PI minscore maxperiod. |
| `discovery.period` | list of int x2 |  | `[80, 300]` | Motif period range [min, max] accepted by Step 3. |
| `discovery.min_copies` | int |  | `10` | Minimum TRF copy number for a motif (Step 3). |
| `discovery.monomerize` | bool |  | `True` | Step 3b: replace a motif that is k tandem copies of a shorter unit (TRF reporting an array at a multiple of its unit, e.g. 292 for a 146-letter unit) by one unit. Test: the motif against itself rotated by len/k, banded global alignment, identity >= monomer_min_id; units never shorter than period[0]. Added 2026-10-18; the published run had no such step (set false to reproduce it exactly). |
| `discovery.monomer_min_id` | float |  | `0.7` | Step 3b: minimum self-identity (motif vs itself rotated by one unit) to call a motif k units. Random and AT-rich random sequence of motif lengths stays far below (0 of 2 000 called in the 2026-10-18 test). |
| `discovery.simple_repeat_cov` | int |  | `50` | Step 4: motifs whose own TRF re-scan covers more than this many bp are dropped as simple repeats. |
| `discovery.strand` | one of `plus`, `both` |  | `plus` | Step 5 canonical rotation and vsearch clustering strand. The published run is plus-only (deliberate; it is what left L541 unclustered). |
| `discovery.cluster_id` | float |  | `0.9` | vsearch identity for per-species clustering (Step 6). |
| `discovery.cross_cluster_id` | float |  | `0.9` | vsearch identity for cross-species clustering (Step 8). |
| `family_assignment` | path | yes |  | Curated TSV: `cluster<TAB>family[<TAB>note]`, cluster = vsearch cluster number from Step 8. Validated at stage 09: every listed cluster must exist, every family must have at least one cluster; unassigned clusters are reported with sizes and best-consensus bits for the curator. |
| `expansion.min_bits` | float |  | `50` | --min-bits: minimum ssearch36 bits to accept a window. |
| `expansion.min_aln_frac` | float |  | `0.67` | --min-aln-frac: fraction of the window (win_mono x MONO) that must align. |
| `expansion.step_mono` | int |  | `3` | --step-monomers: dense-walk step, in monomers. |
| `expansion.win_mono` | int |  | `3` | --win-monomers: probe window width, in monomers. |
| `expansion.query_repeat` | int |  | `2` | --query-repeat N: the ssearch36 query is N+1 tandem copies of the consensus. |
| `expansion.probe_mono` | int |  | `100` | --probe-monomers: coarse-probe initial step, in monomers. |
| `expansion.probe_min_step` | int |  | `50` | --probe-min-step: coarse-probe minimum step, in monomers. |
| `expansion.max_steps` | int |  | `500` | --max-steps: dense-walk cap per direction (0 = unlimited). |
| `expansion.to_scaffold_ends` | bool |  | `False` | --to-scaffold-ends. |
| `expansion.jobs` | int |  | `16` | --jobs: GNU parallel workers. |
| `scan.win_mono` | int |  | `2` | Scan window width in monomers (--scan-win-bp = win_mono x MONO). |
| `scan.step_mono` | int |  | `1` | Scan step in monomers (--scan-step-bp = step_mono x MONO). |
| `scan.min_bits` | float |  | `0` | --scan-min-bits: minimum bits to emit a scan row. |
| `arbitration.delta` | float |  | `20` | DELTA: seed/rescan contradiction threshold (UNIQUE_SEED_CONTRADICTION at delta <= -DELTA). |
| `arbitration.delta_ambig` | float |  | `15` | DELTA_AMBIG: ambiguity flag threshold. |
| `arbitration.seg_bin` | int |  | `1000` | SEG_BIN: D3b bin width in bp. |
| `arbitration.seg_min_run` | int |  | `10` | SEG_MIN_RUN: minimum run of bins to start a segment. |
| `arbitration.seg_min_locus` | int |  | `2000` | SEG_MIN_LOCUS: segments shorter than this are handled per `short_segments`. |
| `arbitration.short_segments` | one of `drop`, `merge` |  | `drop` | SEG_SHORT: `drop` is the published behaviour; `merge` folds short segments into the longer neighbour. |
| `arbitration.bin_min_bits` | float |  |  | BIN_MIN_BITS: a 1-kb bin counts as satellite when its best window reaches this many bits (used by the bin-majority family rule). Absolute scores depend on the monomer length and the scoring, so set it per satellite. Default: expansion.min_bits. |
| `arbitration.bin_share` | float |  | `0.75` | BIN_SHARE: segments, seed conflicts and unseeded loci take the family that wins most 1-kb bins (bins whose best window reaches arbitration.bin_min_bits); the call is _AMBIGUOUS when that family holds less than this share of the bins. Not in the published run, which used the locus-wide maximum window (author's decision 2026-10-25). |
| `orthology.flank_bp` | int |  | `600` | FLANK: flank length extracted outward from each array edge. |
| `orthology.guard_bp` | int |  | `0` | GUARD: extra distance a flank must keep from another locus. |
| `orthology.min_aln` | int |  | `200` | Minimum aligned length for an ortholog hit. |
| `orthology.min_pid` | float |  | `80.0` | Minimum percent identity for an ortholog hit. |
| `subfamily.consensuses` | path |  |  | Directory of per-sub-subfamily consensus FASTA files (`subx/` in the Darevskia run) or one multi-record FASTA. Optional: without it the stage is skipped. |
| `subfamily.members` | map -> list |  |  | Family -> list of its sub-subfamily names (the FAMSUBS table). Each family listed must exist in `families`. |
| `subfamily.confident_delta` | float |  | `10` | Bits margin over the second-best sub-subfamily for a confident call. |
| `terminal.window_mono` | float |  | `1.5` | Window into the array from each terminus, in monomers. |
| `terminal.edge_tol` | int |  | `10` | EDGE_TOL: tolerated gap between alignment and array edge, bp. |
| `finalisation.large_scaffold_bp` | int |  | `50000000` | Scaffolds longer than this are TRF-scanned in overlapping chunks (the published run needed this for one 85.8 Mb unp scaffold). |
| `finalisation.chunk_bp` | int |  | `2000000` | Chunk size for the chunked TRF scan. |
| `finalisation.overlap_bp` | int |  | `50000` | Chunk overlap for the chunked TRF scan. |
| `finalisation.trf_timeout` | int |  | `1800` | Seconds after which a chunk's TRF run is stopped and the chunk bisected (the published run's chunk 44 of the 85 Mb scaffold never finished: a telomeric tail). |
| `finalisation.min_piece_bp` | int |  | `20000` | Pieces shorter than this are no longer bisected; if TRF still times out on them they are listed in unscanned.tsv and left out (never imputed). |
| `finalisation.edge_tol` | int |  | `50` | EDGE_TOL used by edge refinement. |
| `finalisation.rearbitrate` | bool |  | `True` | No effect (kept so old configs load): stage 14 always arbitrates the final catalogue; there is no separate re-arbitration pass. |
| `anchor.family` | str |  |  | Family whose consensus defines the frame. Default: the first family under `families`. |
| `anchor.motif` | str |  | `` | Motif the anchor frame must start with (Darevskia CLsat: AAAGAAG). Empty: the anchor family's consensus is used exactly as given. |
| `anchor.independent` | bool |  | `False` | Families are unrelated (no common frame exists, e.g. different satellites of one genome): each family consensus is its own frame and only sub-subfamilies are rotated, onto their family. Default false: all families must be placeable on the anchor. |
| `anchor.max_motif_mismatch` | int |  | `1` | Mismatches tolerated between the motif and the start of a rotated sequence (diverged families, e.g. CLsat4, may differ by one base). |
| `anchor.min_identity` | float |  | `0.6` | Minimum identity of the local alignment to the anchor for a sequence to be placed; below it the stage fails instead of guessing a frame. |
| `manual_loci` | path |  |  | Optional BED of hand-added loci (L541 in the Darevskia catalogue): chrom start end `<genome>|...` [score strand family]; col 4 must start with a genome code, col 7 (optional) names the family, else the best-scoring family is used. Entries are tagged MANUAL in provenance and listed in the manifest. |
| `tools` | map -> path |  |  | Optional executable overrides: tool name -> path. The special key `bin_dir` is a directory searched before PATH for every tool. |
| `viewer.title` | str |  |  | Page title. |
| `viewer.species_names` | map -> str |  |  | Code -> display name. Keys must exist in `genomes`. |
