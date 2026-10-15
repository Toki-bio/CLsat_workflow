# Viewer

The interactive viewer and its data products live in their own repository and are served from
its root: **[Toki-bio/CLsat-viewer](https://github.com/Toki-bio/CLsat-viewer)** →
<https://toki-bio.github.io/CLsat-viewer/>.

The copy that used to sit here (`index.html` reporting `v2.10-fix-msa-display-20260220` but built on
the 626-locus catalogue, with two empty JSON files and no bitscore tracks) was an unlabelled
intermediate that could not run on its own; it was removed on 2026-10-15. Use the viewer repository
at the commit named in the manuscript's data-availability statement.

Data products the viewer reads, and the scripts here that make them:

| File | Producer |
|---|---|
| `MASTER.merged.bed`, `MASTER.final_assignment.clean.tsv`, `MASTER.locus_bits.tsv` | `03_master_arbitration/` (+ hand-added L541, + `07_post_catalogue/candidate_scan.py` for C1–C7) |
| `MASTER_CLsat{1..4}.bits.tsv.gz`, `MASTER_flank.bits.tsv.gz` | the symmetric 2nd scan (`clisat_expand_then_scan_merged.sh --scan-only`), gzip -9 |
| `flanks.json`, `flank_orthologs.json` | `04_flank_orthology/rebuild_ortho_segmented.py` |
| `subfamily.json`, `terminal_monomers.json` | `05_subfamily_terminal_boundary/assign_subfamily.sh`, `terminal_monomers.sh`, `build_subfam_term_json.py` |
