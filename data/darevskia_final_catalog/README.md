# Final Darevskia CLsat catalog and DarSat tables (deposited data of the CLsat manuscript)

Files of the final run `run_gb_v3` of this workflow (manifest code commit
e45eeb226488cfde1da9ab6c4e4d05be598d7531) on the six Darevskia assemblies val, mix, unv, unn, armf and armm, and
of the DarSat family search on the same panel. They are copied unchanged from the run directories; only the three
DarSat files were renamed. Per-window bitscore tracks, per-locus bitscores and the JSON files of the viewer are in
https://github.com/Toki-bio/CLsat-viewer (`current/`). The single-copy flank alignment sets are in
`../darevskia_single_copy_flank_sets/`.

The final catalog has 594 loci (`locus_master.tsv`). `MASTER.loci.segmented.bed` and `MASTER.final_assignment.tsv`
were written before the haplotype-duplicate rule and still hold L134 (unv, 47,681 bp), an array assembled in both
parental assemblies of D. unisexualis animal 224 and removed by that rule (`hapdup.tsv`, verdict REMOVED).

## CLsat catalog (run_gb_v3)

| file | content | columns |
|---|---|---|
| `MASTER.loci.segmented.bed` | segmented loci (595 lines, incl. L134) | contig, start (0-based), end, locus, assembly |
| `MASTER.final_assignment.tsv` | family arbitration per locus (595 lines, incl. L134) | locus, family, family bits, runner-up family, runner-up bits, margin, arbitration flag, families whose expansions formed the locus |
| `locus_master.tsv` | the final 594 loci | locus, contig, start (0-based), end, assembly, family, monomer length (bp) |
| `locus_subfamily.tsv` | sub-subfamily per locus, called by the whole array in 1-kb pieces (header line) | locus, major, subfamily, pct_id, aln_len, bits, second, second_bits, delta, share (of the length won), counted_bp, n_pieces, runner_up, runner_share |
| `hapdup.tsv` | haplotype-duplicate test of array-only contigs (header line) | locus, genome, contig, contig_len, partner, partner_hit, qcov, ident_gc, partner_loci, verdict |

Arbitration flags: UNIQUE_SEED_WITH_COMPETITOR, UNIQUE_SEED_CONTRADICTION, SEED_CONFLICT_RESOLVED(_AMBIGUOUS),
SEGMENT_RESOLVED(_AMBIGUOUS), NO_SEED_RESOLVED, NO_SEED_AMBIGUOUS (main text Methods 4.4).

## Read check of the catalog (Additional file 1: Table S5)

`read_check_labels.tsv`: one row per locus of the final catalog (594). Each animal's own HiFi reads were mapped to its assembly or assemblies (minimap2 map-hifi; both assemblies of a hybrid together). Columns: locus, assembly, contig, start (0-based), end, family, length, depth_vs_genome_median (read depth over the array / genome median), reads_crossing_left_end, reads_crossing_right_end (mapping quality >= 20), reads_spanning_array, label. Labels (first that applies): proven_whole (a read spans the array, or spanning reads give a length within 20%); not_checkable (both ends at contig ends); flag_end (an end inside the contig crossed by fewer than 2 reads); flag_high / flag_low (depth above 1.6 / below 0.65 times the median); proven_ends (otherwise). An end counts as a contig end only when the array reaches the contig's first or last base. Table S5 gives the share of CLsat length per label and assembly.

## DarSat (darsat_gb)

| file | content |
|---|---|
| `DarSat_consensus_sequences.fasta` | the 39 DarSat family consensuses (= Additional file 4 of the manuscript; md5 ac525483feb173994c24c87e004d8a05) |
| `DarSat_family_by_assembly.tsv` | family, assembly, number of filtered TRF motifs assigned (5,551 in all) |
| `DarSat_family_by_outgroup_genome.tsv` | genome (GenBank accession; CTRL_val = the D. valentini motif set searched with the outgroup command, a control), family, number of motifs assigned |

The search and thresholds are described in Additional file 3 (DarSat Supplementary Methods).
