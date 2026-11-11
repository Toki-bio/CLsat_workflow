#!/bin/bash
set -euo pipefail
# Sub-subfamily assignment from the WHOLE array (subfamily.method: bins; owner's decision 2026-11-11).
# assign_subfamily.sh reads one 2-monomer window at the array boundary; here every locus is cut into 1-kb pieces,
# every piece is searched against the sub-subfamily consensuses of the locus's family, and the locus takes the
# sub-subfamily that wins the most piece length.
#
# Usage: assign_subfamily_bins.sh <work_dir> <subfamily_consensus_dir> <members.tsv> <genomes.tsv> <jobs>
#   work_dir                 holds locus_master.tsv (locus scf start0 end copy family mono);
#                            receives locus_subfamily.tsv and subfamily_bins.tsv
#   subfamily_consensus_dir  one FASTA per sub-subfamily, <name>.fa, DOUBLED (the files assign_subfamily.sh reads)
#   members.tsv              family<TAB>sub-subfamily, one row per sub-subfamily
#   genomes.tsv              copy code<TAB>genome FASTA (indexed, or indexable by samtools faidx)
#   jobs                     ssearch36 threads (-T)
# Tools: samtools, ssearch36, gawk on PATH.
#
# Rule:
#   pieces   piece i of a locus = [i, i+1000) from the locus start, i = 0, 1000, 2000, ...; the pieces stop at the
#            first piece shorter than 300 bp when i > 0 (a short last piece is dropped; a locus shorter than 1 kb is
#            one piece).
#   search   ssearch36 -n -z 3 -m 8C -E 1000 -T <jobs> <pieces of family F> <library of F>; the library holds every
#            record of <name>.fa for each sub-subfamily of F in members.tsv order, named by the sub-subfamily.
#   piece    best = the sub-subfamily with the highest bits over all hits of the piece (the maximum per sub-subfamily
#            is kept); second = the next one. Equal bits: the sub-subfamily that ssearch36 reported first for the
#            piece. A piece with no hit counts for nobody.
#   locus    call = the sub-subfamily with the largest total length of the pieces it won; runner-up = the second;
#            share = length won / length of the pieces that had a hit.
#            TIE RULE: equal lengths -> the sub-subfamily whose first won piece lies nearest the locus start (the
#            order in which the sub-subfamilies are first reached along the locus). Ties are counted in the log.
#   A family with one sub-subfamily gets that call wherever a piece has a hit. A family without members.tsv rows,
#   and a locus none of whose pieces has a hit, get NA.
#
# locus_subfamily.tsv (one row per locus, locus_master.tsv order). The first 9 columns keep the names and positions
# of assign_subfamily.sh so that its readers keep working, with these meanings:
#   locus major subfamily(call) pct_id(NA) aln_len(counted bp) bits(share x 100, 1 decimal)
#   second(runner-up or NA) second_bits(runner-up share x 100) delta(bits - second_bits, percentage points; NA
#   without a runner-up)   then: share counted_bp n_pieces runner_up runner_share
# subfamily_bins.tsv (one row per piece): copy locus start end family best bits second second_bits
#   (start/end = offsets inside the locus, 0-based half-open; NA = no hit / no second sub-subfamily).
[[ $# -eq 5 ]] || { echo "Usage: $0 <work_dir> <subfamily_consensus_dir> <members.tsv> <genomes.tsv> <jobs>" >&2; exit 1; }
WORK="$1"; SUBX="$2"; MEMBERS="$3"; GENOMES="$4"; JOBS="$5"
PIECE_BP=1000          # piece length
MIN_LAST_BP=300        # a last piece shorter than this is dropped
MASTER=$WORK/locus_master.tsv
GAWK=${GAWK:-gawk}
command -v samtools >/dev/null || { echo "samtools not in PATH" >&2; exit 1; }
command -v ssearch36 >/dev/null || { echo "ssearch36 not in PATH" >&2; exit 1; }
command -v "$GAWK" >/dev/null || { echo "gawk not in PATH" >&2; exit 1; }
[[ -d "$SUBX" ]] || { echo "sub-subfamily consensus dir not found: $SUBX" >&2; exit 1; }
[[ -s "$MASTER" ]] || { echo "missing $MASTER" >&2; exit 1; }
[[ -f "$MEMBERS" ]] || { echo "missing $MEMBERS" >&2; exit 1; }
[[ -f "$GENOMES" ]] || { echo "missing $GENOMES" >&2; exit 1; }
[[ "$JOBS" =~ ^[1-9][0-9]*$ ]] || { echo "jobs must be a positive integer: $JOBS" >&2; exit 1; }
export LC_ALL=C

cd "$WORK"
TMPD=$(mktemp -d "$WORK/subfam_bins.XXXXXX")
trap 'rm -rf -- "$TMPD"' EXIT

# ---- 1. locus sequences, one samtools call per genome copy -> seq.tsv (locus<TAB>sequence)
declare -A GENOME
while IFS=$'\t' read -r code path _; do
  [[ -z "$code" || "$code" == \#* ]] && continue
  GENOME[$code]=$path
done < "$GENOMES"
: > "$TMPD/seq.tsv"
n=0
while IFS= read -r sp; do
  n=$((n + 1))
  [[ -n "${GENOME[$sp]:-}" ]] || { echo "FATAL: copy $sp of locus_master.tsv is not in $GENOMES" >&2; exit 1; }
  "$GAWK" -F'\t' -v S="$sp" '$5 == S { print $2 ":" $3 + 1 "-" $4 }' "$MASTER" > "$TMPD/reg.$n"
  samtools faidx "${GENOME[$sp]}" -r "$TMPD/reg.$n" > "$TMPD/reg.$n.fa" \
    || { echo "FATAL: samtools faidx failed for copy $sp (${GENOME[$sp]})" >&2; exit 1; }
  # records come back in the order asked; a wrong count or length is an error, never a silent shift
  "$GAWK" -F'\t' -v S="$sp" -v FA="$TMPD/reg.$n.fa" '
    BEGIN { r = 0; while ((getline l < FA) > 0) { if (l ~ /^>/) r++; else seq[r] = seq[r] l } close(FA) }
    $5 == S { k++; if (length(seq[k]) != $4 - $3) { printf "FATAL: %s (%s:%d-%d): got %d letters, expected %d\n", $1, $2, $3 + 1, $4, length(seq[k]), $4 - $3 > "/dev/stderr"; bad = 1; exit 1 }
              print $1 "\t" seq[k] }
    END { if (!bad && k != r) { printf "FATAL: copy %s: %d loci, %d sequences returned\n", S, k, r > "/dev/stderr"; exit 1 } }
  ' "$MASTER" >> "$TMPD/seq.tsv"
done < <("$GAWK" -F'\t' '!($5 in s) { s[$5] = 1; print $5 }' "$MASTER")

# ---- 2. pieces: pieces.tsv (copy locus start end family, locus_master.tsv order) and one query FASTA per family
"$GAWK" -F'\t' -v P="$PIECE_BP" -v M="$MIN_LAST_BP" -v D="$TMPD" '
  FNR == NR { seq[$1] = $2; next }
  {
    if (!($6 in fn)) { fn[$6] = ++nf; print $6 > (D "/families") }
    s = seq[$1]; L = length(s); q = D "/q." fn[$6] ".fa"
    for (i = 0; i < L; i += P) {
      b = substr(s, i + 1, P); lb = length(b)
      if (lb < M && i > 0) break
      printf ">%s|%s|%d|%d\n%s\n", $5, $1, i, i + lb, b > q
      printf "%s\t%s\t%d\t%d\t%s\n", $5, $1, i, i + lb, $6 > (D "/pieces.tsv")
    }
  }' "$TMPD/seq.tsv" "$MASTER"
: >> "$TMPD/pieces.tsv"; : >> "$TMPD/families"

# ---- 3. one search per family
: > "$TMPD/hits.m8"
n=0
while IFS= read -r fam; do
  n=$((n + 1))
  [[ -s "$TMPD/q.$n.fa" ]] || continue
  lib=$TMPD/lib.$n.fa; : > "$lib"
  while IFS= read -r sub; do
    consx=$SUBX/${sub}.fa
    # every listed sub-subfamily must have its consensus: a missing one would silently narrow the choice
    [[ -s "$consx" ]] || { echo "FATAL: sub-subfamily $sub of $fam (members.tsv) has no consensus file $consx" >&2; exit 1; }
    "$GAWK" -v N="$sub" '
      { sub(/\r$/, "") }
      /^>/ { if (h) printf ">%s\n%s\n", N, s; h = 1; s = ""; next }
      h { s = s $0 }
      END { if (h) printf ">%s\n%s\n", N, s }' "$consx" >> "$lib"
  done < <("$GAWK" -F'\t' -v F="$fam" '$1 == F && $2 != "" { print $2 }' "$MEMBERS")
  if [[ ! -s "$lib" ]]; then echo "WARNING: family $fam has no sub-subfamily in $MEMBERS: its loci get NA" >&2; continue; fi
  ssearch36 -n -z 3 -m 8C -E 1000 -T "$JOBS" "$TMPD/q.$n.fa" "$lib" > "$TMPD/hits.$n.m8" \
    || { echo "FATAL: ssearch36 failed for family $fam" >&2; exit 1; }
  cat "$TMPD/hits.$n.m8" >> "$TMPD/hits.m8"
done < "$TMPD/families"

# ---- 4. best sub-subfamily per piece, call per locus
OUT=$WORK/locus_subfamily.tsv
BINS=$WORK/subfamily_bins.tsv
"$GAWK" -F'\t' -v OUT="$OUT" -v BINS="$BINS" -v TIES="$TMPD/ties" '
  # pick(n, val, skip): index of the largest val among 1..n (the first one on equal values), skipping index skip
  function pick(n, val, skip,   i, b) { b = 0; for (i = 1; i <= n; i++) if (i != skip && (b == 0 || val[i] > val[b])) b = i; return b }
  BEGIN {
    printf "locus\tmajor\tsubfamily\tpct_id\taln_len\tbits\tsecond\tsecond_bits\tdelta\tshare\tcounted_bp\tn_pieces\trunner_up\trunner_share\n" > OUT
    printf "copy\tlocus\tstart\tend\tfamily\tbest\tbits\tsecond\tsecond_bits\n" > BINS
    ptie = 0; ltie = 0
  }
  FILENAME == ARGV[1] {                                   # ssearch36 -m 8C rows: query subject ... bits (column 12)
    if ($0 == "" || $0 ~ /^#/) next
    b = $12 + 0
    if (!(($1, $2) in hb)) { hn[$1]++; ho[$1, hn[$1]] = $2; hb[$1, $2] = b; ht[$1, $2] = $12 }
    else if (b > hb[$1, $2]) { hb[$1, $2] = b; ht[$1, $2] = $12 }
    next
  }
  FILENAME == ARGV[2] {                                   # pieces, in locus order
    q = $1 "|" $2 "|" $3 "|" $4; n = hn[q] + 0
    delete v; for (i = 1; i <= n; i++) v[i] = hb[q, ho[q, i]]
    b1 = pick(n, v, 0); b2 = b1 ? pick(n, v, b1) : 0
    s1 = b1 ? ho[q, b1] : "NA"; s2 = b2 ? ho[q, b2] : "NA"
    t1 = b1 ? ht[q, s1] : "NA"; t2 = b2 ? ht[q, s2] : "NA"
    gsub(/ /, "", t1); gsub(/ /, "", t2)
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, $5, s1, t1, s2, t2 > BINS
    np[$2]++
    if (b1) {
      if (b2 && v[b1] == v[b2]) ptie++
      if (!(($2, s1) in won)) { wn[$2]++; wo[$2, wn[$2]] = s1 }
      won[$2, s1] += $4 - $3
    }
    next
  }
  {                                                       # locus_master.tsv
    L = $1; n = wn[L] + 0; tot = 0
    delete v; for (i = 1; i <= n; i++) { v[i] = won[L, wo[L, i]]; tot += v[i] }
    c1 = pick(n, v, 0); c2 = c1 ? pick(n, v, c1) : 0
    if (!c1) { printf "%s\t%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t0.000\t0\t%d\tNA\t0.000\n", L, $6, np[L] + 0 > OUT; next }
    p1 = sprintf("%.1f", 100 * v[c1] / tot)
    if (c2) {
      p2 = sprintf("%.1f", 100 * v[c2] / tot)
      if (v[c1] == v[c2]) { ltie++; printf "%s\t%s\t%s\t%d\n", L, wo[L, c1], wo[L, c2], v[c1] > TIES }
      printf "%s\t%s\t%s\tNA\t%d\t%s\t%s\t%s\t%.1f\t%.3f\t%d\t%d\t%s\t%.3f\n", L, $6, wo[L, c1], tot, p1, wo[L, c2], p2, p1 - p2, v[c1] / tot, tot, np[L], wo[L, c2], v[c2] / tot > OUT
    } else
      printf "%s\t%s\t%s\tNA\t%d\t%s\tNA\tNA\tNA\t%.3f\t%d\t%d\tNA\t0.000\n", L, $6, wo[L, c1], tot, p1, v[c1] / tot, tot, np[L] > OUT
  }
  END {
    printf "pieces with equal best and second bits (first reported wins): %d\n", ptie
    printf "loci with equal won length for call and runner-up (tie rule applied): %d\n", ltie
  }' "$TMPD/hits.m8" "$TMPD/pieces.tsv" "$MASTER"

echo "rows: $(( $(wc -l < "$OUT") - 1 ))   pieces: $(( $(wc -l < "$BINS") - 1 ))"
if [[ -s "$TMPD/ties" ]]; then echo "=== tied loci (locus call runner_up won_bp) ==="; cat "$TMPD/ties"; fi
echo "=== sample ==="
head -8 "$OUT"
echo "=== subfamily distribution ==="
"$GAWK" -F'\t' 'NR>1{print $3}' "$OUT" | sort | uniq -c | sort -rn
