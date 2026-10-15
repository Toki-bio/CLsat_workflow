# legacy — superseded code, kept for provenance only

Nothing here is current. It is retained because the manuscript's history refers to it.

- `superseded_blastn/`: the earlier blastn-based flank-orthology method and its consumers
  (`disc_ortho.py`, `flank_classify.py`, `ortho_matrix.py`). It lacks the flank overlap guard and
  compares same-side flanks only, giving 45 pairs instead of the correct 53. The authoritative
  method is `satpipe/stages/04_flank_orthology/clean_ortho_pairs.py` /
  `rebuild_ortho_segmented.py`. See `superseded_blastn/README.md`.
- `expand_clsats.sh`: the serial predecessor of `clisat_expand_then_scan_merged.sh`, with the
  `--min-aln-frac` base bug (default 1.0 against the query length made the pass condition
  unsatisfiable and truncated CLsat1 walks).
