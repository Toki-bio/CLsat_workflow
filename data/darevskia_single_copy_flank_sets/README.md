# CLsat single-copy flank sets (deposited data of the CLsat manuscript)

One set per single-copy flank of the final CLsat catalogue: 53 flanks, one each on 53 of the 594 loci
(no locus has two). These are the flanks that the workflow uses for cross-assembly orthology
(main text, Methods section 5; Additional file 2, Step 13c).

## What the query is
Each query is the flank window exactly as the single-copy test saw it: 600, 3,000 or 5,000 bp outward from
the array edge, with the parts that the workflow masked written as N (pieces matching a CLsat consensus,
tandem repeats, and pieces with another copy in the own genome). Only the unmasked letters were tested and
found single-copy. The query header gives copy, contig, window coordinates (1-based), masked and unmasked bp.

## How the matches were found
The query was searched against each of the six genome copies (val, mix, unv, unn, armf, armm) with the
minimap2 command of the workflow's single-copy test (minimap2 2.28; -k13 -w4 -s40 -m25 -A2 -B4 -O4,24 -E2,1
-c -N 200 -p 0). A match is kept when it has at least 200 aligned letters and at least 80% identity
(matching letters / alignment block), the thresholds of the test. Up to 60 matches per flank are given,
ordered by matching letters.

## Files
- `<locus>_<L|R>.fa`      query first, then the matched genomic pieces in the query's orientation.
- `<locus>_<L|R>.aln.fa`  MAFFT (--adjustdirection --auto) alignment of the same.
- `manifest.tsv`          one row per flank (columns below).

## Match headers
`<copy>_<contig>_<start>-<end>_<plus|minus>_pid<NN>_<tag>`, 1-based coordinates. Tags:
- `SELF` the flank's own place;
- `NEARLOCUS_<locus>_<bp>` within one window length of a catalogue CLsat locus of that copy
  (candidate orthologous position);
- `nearest<locus>_<bp>` farther from any locus; `nolocus` on a contig without a catalogue locus.

## manifest.tsv columns
flank, locus, side, src_copy, window_bp, unmasked_bp, total_hits, self_hit (alignment blocks on the
flank's own place; 2 when the own place aligns in two pieces), other_own_copy
(matches elsewhere in the own copy; 0 for a single-copy flank), other_copies (matches per other copy),
shown, n_near_locus.
