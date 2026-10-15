# Files produced by each stage

The stage contracts from HANDOFF.md, as files. Coordinate frame is stated for every file that
carries positions. Columns marked *to document* are read from the scripts in work-plan step 5
when the driver is wired; until then this table is the index, not yet the full dictionary.

Frames: **1-based inclusive** = TRF `.dat`, Step-3 headers, `samtools faidx` regions, ssearch36
`-m 8C` columns 7-10, scan `pos`. **0-based half-open** = every BED, every locus `start`/`end` in
JSON. Conversions live in `make_bed_from_clsats_headers.sh:25`, `clisat_expand_then_scan_merged.sh`
(263-297, 559-579) and `arbitration_patched.sh` D3/D3b (`pos-1`).

| stage | file | format | frame | notes |
|---|---|---|---|---|
| 02 | `<code>.dat` | TRF `.dat` (`Sequence:` header lines, then `start end period copies ... motif sequence`) | 1-based | `trf <genome> 2 5 7 80 10 50 2000 -h`; exit status is not 0 on success in every build |
| 03 | `motifs_output/<code>.fa` | FASTA; header `code|SCF:START-END|len=P|copies=C|locuslen=N` | 1-based | one record per TRF motif with period in `discovery.period` and copies >= `min_copies`; sequence uppercased |
| 04 | `<code>.fa.coverage.txt` | `header<SPACE>covered_bp` | | internal TRF coverage of each motif |
| 04 | `<code>.fa.clean.fas`, `<code>.fa.repeats.fas` | FASTA | | split at `simple_repeat_cov` (<= 50 bp clean) |
| 05 | `clean/<code>.fa.clean.fas.rot` | FASTA | | lexicographically minimal rotation, plus strand only |
| 06 | `clusters_by_file/<code>/<code>_cluster_N.fa` | FASTA | | vsearch `--cluster_fast --id 0.90` (strand plus); N = vsearch cluster number + 1 |
| 07 | `<cluster>_rotated/<cluster>.aligned.fasta`, `<cluster>.rotated.fa` | aligned FASTA (gaps) | | MAFFT L-INS-i; rotated to the rightmost first non-gap column, realigned |
| 08 | `all7_pooled.degap.fa`, `cross_species_centroids.fa`, `cross_species.uc`, `cross_clusters/clust_N` | FASTA / vsearch UC | | references prepended so they win centroid election; N is the cluster number used in `family_assignment` |
| 09 | `<family>.members.fa` | FASTA | 1-based (headers) | genomic members of the family's clusters; reference records dropped |
| 10 | `<family>.bed` | BED6 `chrom start end name 0 .` | 0-based | from the `SCF:START-END` in each header; name = full header |
| 11 | `out_<family>.expanded_per_locus.bed` | BED4 | 0-based | one expanded interval per seed (empty with `--scan-only`) |
| 11 | `out_<family>.expanded_per_locus.tsv` | `species chrom start0 end0 name seed_start0 seed_end0 reason_L reason_R expand_L_bp expand_R_bp expand_total_bp` | 0-based | stop reasons: `signal_lost flip past_scaffold_start past_scaffold_end window_invalid max_steps coarse_probe` |
| 11 | `out_<family>.expanded_merged.bed` | BED4 | 0-based | merged per species+scaffold; input to stage 12 |
| 11 / 13 | `*.bits.tsv` | `scaffold pos bits name species` | `pos` 1-based | scan profile; use `pos-1` as BED start |
| 11 | `out_<family>.dropped.tsv` | *to document* | | every seed that could not be processed, with a reason; seed count reconciled |
| 12 | `MASTER.merged.bed` | BED5 `chrom start end tag provenance` | 0-based | tag `code|SCF:START1-END|MASTER` (START1 is 1-based); provenance = comma-joined families |
| 13 | `MASTER_<family>.bits.tsv` | as `*.bits.tsv` | `pos` 1-based | symmetric rescan, window 2xMONO, step MONO |
| 14 | `MASTER.loci.bed` | BED4, name `L<n>` | 0-based | `L"NR` over the `LC_ALL=C`-sorted MASTER |
| 14 | `MASTER.locus_provenance.tsv` | `locus provenance` | | |
| 14 | `MASTER.loci.segmented.bed` (+ `.short_segments.bed`) | BED4, names `L<n>` or `L<n>.<k>` | 0-based | D3b segmentation; short segments dropped (`short_segments: drop`) or merged |
| 14 | `MASTER.locus_seed_hits.tsv`, `MASTER.locus_seeds.tsv` | *to document* | | seed overlap per locus; `SEED_CONFLICT` marked |
| 14 | `MASTER.<family>.maxbits.tsv`, `MASTER.locus_bits.tsv` | `locus bits_fam1 ... bits_famN` | | per-locus maxima of the rescan |
| 14 | `MASTER.final_assignment.tsv` | `locus_id family assigned_bits best_other best_other_bits delta flag provenance` (8 columns, no header in the `.clean.tsv` the viewer reads) | | delta = assigned - strongest other for every flag; flags: `UNIQUE_SEED_ONLY UNIQUE_SEED_WITH_COMPETITOR UNIQUE_SEED_CONTRADICTION SEED_CONFLICT_RESOLVED[_AMBIGUOUS] SEGMENT_RESOLVED[_AMBIGUOUS] NO_SEED_RESOLVED NO_SEED_AMBIGUOUS`; viewer adds `CANDIDATE_UNCATALOGUED` |
| orthology | `flanks.json` | JSON, per locus | 1-based for L-loci (from `rebuild_ortho_segmented.py`), 0-based for C-loci (from `candidate_scan.py`) | fix the C-locus frame when touched |
| orthology | `flank_orthologs.json`, matrix | JSON / TSV | | pairs counted per locus pair (53), not per flank link (76); best hit per target assembly applies to the matrix |
| subfamily | `manifest.tsv` | `locus kind sp scf rs re trunc explen` (kinds LF/LB/RB/RF) | 1-based | from `run_extract.sh`; `all_regions.fa` holds the sequences |
| subfamily | `locus_subfamily.tsv` | `locus major subfamily pct_id aln_len bits second second_bits delta` | | bitscores compared as floats |
| subfamily | `subfamily.json` | JSON, per locus `{major, subfamily, pid, bits, second, delta, confident}` | | confident if delta >= `confident_delta` or no competitor |
| subfamily | `terminal_monomers.tsv` / `.json` | `cons_start cons_end covered partial strand edge_gap pct_id` per terminus | consensus coordinates 1-based | rewritten 2026-10-15; phase read at the array edge |
| finalisation | `chunk_coords.tsv`, `chunks/chunk_N.{fa,dat}` | | 1-based | 2 Mb chunks, 50 kb overlap, one scaffold |
| finalisation | `refined_*.bed`, `.report.tsv` | *to document* | 0-based | nucleotide-level edge refinement |
| viewer | `scaffold_lengths.tsv`, `<family>.cons.fa`, `MASTER_*.bits.tsv.gz`, `MASTER_flank.bits.tsv.gz` | | | inputs the viewer reads; gzip -9 -n |
