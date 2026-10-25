#!/bin/bash
set -euo pipefail
# Build locus_master.tsv (Lid scf start0 end species family MONO) from the stage-12/14 outputs.
#
# Usage: build_locus_master.sh <MASTER.merged.bed> <arbitration_out_dir> <work_dir> <families.tsv>
#   MASTER.merged.bed    not read (kept for the calling convention); the genome comes from col 5 of
#                        MASTER.loci.segmented.bed
#   arbitration_out_dir  holds MASTER.final_assignment.tsv and MASTER.loci.segmented.bed (stage 14)
#   work_dir             created; receives locus_master.tsv (+ helper tables)
#   families.tsv         family  mono  [consensus.fa]
[[ $# -eq 4 ]] || { echo "Usage: $0 <MASTER.merged.bed> <arbitration_out_dir> <work_dir> <families.tsv>" >&2; exit 1; }
MASTER_BED="$1"; OUT="$2"; WORK="$3"; FAMTAB="$4"
mkdir -p "$WORK"

# 1. the genome of every locus is col 5 of MASTER.loci.segmented.bed (written by stage 14/15). It is never
#    looked up by contig name: two genomes may share contig names (hifiasm h1tg/h2tg in two animals).
awk -F'\t' 'NF<5 || $5==""{ print "FATAL: MASTER.loci.segmented.bed row without genome (col 5): " $0 > "/dev/stderr"; bad=1; exit }
  END{ exit bad }' "$OUT/MASTER.loci.segmented.bed"

# 2. locus -> family
awk -F'\t' '{print $1"\t"$2}' "$OUT/MASTER.final_assignment.tsv" > "$WORK/locus2subtype.tsv"

# 3. MONO table from families.tsv
grep -v '^#' "$FAMTAB" | awk -F'\t' 'NF>=2{print $1"\t"$2}' > "$WORK/mono.tsv"

# 4a. segmented bed + species
awk -F'\t' '{print $4"\t"$1"\t"$2"\t"$3"\t"$5}' "$OUT/MASTER.loci.segmented.bed" > "$WORK/tmp1.tsv"

# 4b. + family  (use subt[] - 'sub' is a reserved awk function)
awk -F'\t' 'FNR==NR{subt[$1]=$2;next}{print $0"\t"(($1 in subt)?subt[$1]:"NA")}' \
  "$WORK/locus2subtype.tsv" "$WORK/tmp1.tsv" > "$WORK/tmp2.tsv"

# 4c. + MONO
awk -F'\t' 'FNR==NR{mono[$1]=$2;next}{print $0"\t"(($6 in mono)?mono[$6]:0)}' \
  "$WORK/mono.tsv" "$WORK/tmp2.tsv" > "$WORK/locus_master.tsv"
rm -f "$WORK/tmp1.tsv" "$WORK/tmp2.tsv"

echo "=== locus_master.tsv (first 5) ==="
head -5 "$WORK/locus_master.tsv"
echo "total: $(wc -l < "$WORK/locus_master.tsv")"
echo "NA/zero rows: $(awk -F'\t' '$5=="NA"||$6=="NA"||$7==0' "$WORK/locus_master.tsv" | wc -l)"
echo "cols in row 1: $(head -1 "$WORK/locus_master.tsv" | awk -F'\t' '{print NF}')"
echo "family distribution:"
awk -F'\t' '{print $6}' "$WORK/locus_master.tsv" | sort | uniq -c
