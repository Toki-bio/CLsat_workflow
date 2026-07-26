# Superseded — earlier blastn-based flank orthology

**These scripts did not produce the results reported in the manuscript.** They are kept
for provenance only. Do not run them expecting the published numbers.

| Script | Why it was replaced |
|---|---|
| `extract_flanks.sh` | No overlap guard: a 600 bp flank window was allowed to fall inside a neighbouring annotated CLsat locus, so satellite sequence was admitted as "unique" flank. This contaminated roughly 35% of the flanks in the earlier analysis. |
| `flank_ortho_blast.sh` | Used `blastn`, and compared flanks **within each side only** (left-vs-left, right-vs-right). Flanks are oriented array→outward, so two orthologous arrays in opposite orientations legitimately match left-flank-to-right-flank; the same-side restriction silently discards those. It yields **45 pairs** instead of the correct **53**. |
| `build_ortho_json.py` | Consumes the `.blast` tabular output of the script above. |

The current method is in the parent directory:

- **`clean_ortho_pairs.py`** — self-contained and authoritative. Rebuilds loci, extracts
  flanks *with* the overlap guard, classifies them, then compares them all-versus-all with
  `ssearch36` (**not** blastn, and with **no side restriction**), keeping the single best
  cross-assembly hit per target assembly at ≥200 bp aligned and ≥80% identity. Produces
  `clean_flank_class.tsv` and `clean_ortho_matrix.txt` in the data deposit, plus the
  per-pair FASTA alignments distributed as Additional file 4.
- **`rebuild_ortho_segmented.py`** — the same method keyed to the segmented catalogue
  rather than the merged one, producing `segmented627_*` in the data deposit. It returns
  the identical 53-pair / 18-both-flank result, confirming the outcome does not depend on
  which catalogue the flanks are attributed to.
