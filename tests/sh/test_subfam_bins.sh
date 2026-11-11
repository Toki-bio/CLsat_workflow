#!/usr/bin/env bash
# subfamily.method: bins - satpipe/stages/05_subfamily_terminal_boundary/assign_subfamily_bins.sh on a toy catalogue:
# piece cutting (1-kb pieces, a last piece under 300 bp dropped, a locus under 1 kb = one piece), best/second per piece,
# the call per locus, both tie rules, no-hit pieces and loci, a one-member family, a family without members; then the
# reader build_subfam_term_json.py on the 14-column table, and the method switch _p8_subfamily_method (subx.sh).
# No external tools: samtools and ssearch36 are STUBS written here (the real search is checked on the server run).
#   stub ssearch36: bits of a piece against a library record = (count of the record's first letter in the piece) / 10,
#   one row per record with a count > 0, in library order, plus a weaker second row (the other strand).
# Usage: tests/sh/test_subfam_bins.sh
set -euo pipefail
export LC_ALL=C
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
ROOT=$(cd -- "$HERE/../.." && pwd)
SCRIPT=$ROOT/satpipe/stages/05_subfamily_terminal_boundary/assign_subfamily_bins.sh
TMP=$(mktemp -d "${TMPDIR:-/tmp}/satpipe_bins.XXXXXX")
trap 'rm -rf -- "$TMP"' EXIT
N=0
check() { local label=$1; shift; N=$((N + 1)); if "$@"; then printf 'ok   %3d %s\n' "$N" "$label"; else printf 'FAIL %3d %s\n' "$N" "$label"; exit 1; fi; }
same() { cmp -s "$1" "$2" || { diff "$1" "$2" | head -20; return 1; }; }
rep() { gawk -v s="$1" -v n="$2" 'BEGIN { for (i = 0; i < n; i++) printf "%s", s }'; }

# ---- stubs
mkdir -p "$TMP/bin"
cat > "$TMP/bin/samtools" <<'EOF'
#!/usr/bin/env bash
# stub: samtools faidx GENOME -r REGIONS   (single-line FASTA records; regions name:from-to, 1-based inclusive)
[ "$1" = faidx ] && [ "$3" = -r ] || { echo "stub samtools: unexpected arguments: $*" >&2; exit 2; }
gawk 'FNR == NR { if ($0 ~ /^>/) n = substr($0, 2); else s[n] = $0; next }
      { r = $0; c = r; sub(/:[0-9]+-[0-9]+$/, "", c); split(substr(r, length(c) + 2), p, "-")
        if (!(c in s)) { print "stub samtools: no " c > "/dev/stderr"; exit 1 }
        print ">" r; print substr(s[c], p[1], p[2] - p[1] + 1) }' "$2" "$4"
EOF
cat > "$TMP/bin/ssearch36" <<'EOF'
#!/usr/bin/env bash
# stub: ssearch36 -n -z 3 -m 8C -E 1000 -T JOBS QUERY LIB
[ "$*" = "-n -z 3 -m 8C -E 1000 -T 2 ${*: -2:1} ${*: -1}" ] || { echo "stub ssearch36: unexpected arguments: $*" >&2; exit 2; }
gawk 'FNR == NR { if ($0 ~ /^>/) { n++; name[n] = substr($0, 2) } else lib[n] = substr($0, 1, 1); next }
      /^>/ { q = substr($0, 2); print "# SSEARCH stub query " q; next }
      { for (i = 1; i <= n; i++) { c = gsub(lib[i], "&")
          if (c > 0) { printf "%s\t%s\t100.00\t%d\t0\t0\t1\t%d\t1\t%d\t1e-9\t%.1f\n", q, name[i], c, c, c, c / 10
                       printf "%s\t%s\t90.00\t%d\t0\t0\t1\t%d\t%d\t1\t1e-3\t%.1f\n", q, name[i], c, c, c, c / 20 } } }' "${*: -1}" "${*: -2:1}"
EOF
chmod +x "$TMP/bin/samtools" "$TMP/bin/ssearch36"
export PATH=$TMP/bin:$PATH

# ---- toy catalogue
W=$TMP/reg; mkdir -p "$W" "$TMP/subx"
{ printf '>s1\n'; rep T 10
  rep A 1000; rep C 600; rep A 400; rep C 250; rep G 250         # L1  2500 bp: A1 | A2 (A1 second) | A2 = A3 (equal bits)
  rep T 5
  rep G 1000; rep A 1000; rep C 200                               # L2  2200 bp: A3 | A1 | 200 bp dropped
  printf '\n'; } > "$TMP/ca.fa"
{ printf '>s2\n'
  rep A 700                                                       # L3  700 bp, family FB (one member)
  rep T 1000; rep C 300                                           # L4  1300 bp: no hit | A2 (300 bp is kept)
  rep T 400                                                       # L5  no hit at all
  rep A 500                                                       # L6  family FC: no members
  printf '\n'; } > "$TMP/cb.fa"
printf 'ca\t%s\ncb\t%s\n' "$TMP/ca.fa" "$TMP/cb.fa" > "$TMP/genomes.tsv"
printf 'FA\tA1\nFA\tA2\nFA\tA3\nFB\tB1\n' > "$TMP/members.tsv"
printf '>A1\nAAAA\n' > "$TMP/subx/A1.fa"; printf '>A2 x\nCC\nCC\n' > "$TMP/subx/A2.fa"
printf '>whatever\nGGGG\n' > "$TMP/subx/A3.fa"; printf '>B1\nAAAA\n' > "$TMP/subx/B1.fa"
printf '%s\n' $'L1\ts1\t10\t2510\tca\tFA\t146' $'L2\ts1\t2515\t4715\tca\tFA\t146' $'L3\ts2\t0\t700\tcb\tFB\t100' \
              $'L4\ts2\t700\t2000\tcb\tFA\t146' $'L5\ts2\t2000\t2400\tcb\tFA\t146' $'L6\ts2\t2400\t2900\tcb\tFC\t50' > "$W/locus_master.tsv"

bash "$SCRIPT" "$W" "$TMP/subx" "$TMP/members.tsv" "$TMP/genomes.tsv" 2 > "$TMP/log" 2> "$TMP/err" || { cat "$TMP/log" "$TMP/err"; echo "FAIL the script exits 0"; exit 1; }

cat > "$TMP/exp_loc" <<'EOF'
locus	major	subfamily	pct_id	aln_len	bits	second	second_bits	delta	share	counted_bp	n_pieces	runner_up	runner_share
L1	FA	A2	NA	2500	60.0	A1	40.0	20.0	0.600	2500	3	A1	0.400
L2	FA	A3	NA	2000	50.0	A1	50.0	0.0	0.500	2000	2	A1	0.500
L3	FB	B1	NA	700	100.0	NA	NA	NA	1.000	700	1	NA	0.000
L4	FA	A2	NA	300	100.0	NA	NA	NA	1.000	300	2	NA	0.000
L5	FA	NA	NA	NA	NA	NA	NA	NA	0.000	0	1	NA	0.000
L6	FC	NA	NA	NA	NA	NA	NA	NA	0.000	0	1	NA	0.000
EOF
cat > "$TMP/exp_bins" <<'EOF'
copy	locus	start	end	family	best	bits	second	second_bits
ca	L1	0	1000	FA	A1	100.0	NA	NA
ca	L1	1000	2000	FA	A2	60.0	A1	40.0
ca	L1	2000	2500	FA	A2	25.0	A3	25.0
ca	L2	0	1000	FA	A3	100.0	NA	NA
ca	L2	1000	2000	FA	A1	100.0	NA	NA
cb	L3	0	700	FB	B1	70.0	NA	NA
cb	L4	0	1000	FA	NA	NA	NA	NA
cb	L4	1000	1300	FA	A2	30.0	NA	NA
cb	L5	0	400	FA	NA	NA	NA	NA
cb	L6	0	500	FC	NA	NA	NA	NA
EOF
check "locus_subfamily.tsv: calls, shares, tie by first reached, NA rows"  same "$W/locus_subfamily.tsv" "$TMP/exp_loc"
check "subfamily_bins.tsv: pieces, best/second, equal bits -> first reported" same "$W/subfamily_bins.tsv" "$TMP/exp_bins"
check "log: one piece tie"                         grep -q '^pieces with equal best and second bits (first reported wins): 1$' "$TMP/log"
check "log: one locus tie, listed"                 bash -c "grep -q 'tie rule applied): 1\$' '$TMP/log' && grep -q \$'^L2\tA3\tA1\t1000\$' '$TMP/log'"
check "family without members: warning"            grep -q 'family FC has no sub-subfamily' "$TMP/err"
check "no temporary directory left"                bash -c "! ls -d '$W'/subfam_bins.* >/dev/null 2>&1"

# a listed sub-subfamily without its consensus file is an error
printf 'FA\tA1\nFA\tA9\n' > "$TMP/members_bad.tsv"
check "missing consensus file: exit 1"             bash -c "! bash '$SCRIPT' '$W' '$TMP/subx' '$TMP/members_bad.tsv' '$TMP/genomes.tsv' 2 > /dev/null 2> '$TMP/err2' && grep -q 'A9 of FA' '$TMP/err2'"
# a locus beyond the scaffold end (fewer letters than asked) is an error
W2=$TMP/reg2; mkdir -p "$W2"; printf 'L1\ts2\t2400\t3000\tcb\tFB\t100\n' > "$W2/locus_master.tsv"
check "short sequence: exit 1"                     bash -c "! bash '$SCRIPT' '$W2' '$TMP/subx' '$TMP/members.tsv' '$TMP/genomes.tsv' 2 > /dev/null 2> '$TMP/err3' && grep -q 'got 500 letters, expected 600' '$TMP/err3'"

# ---- the reader of locus_subfamily.tsv takes the 14-column table (first 9 columns)
PY=${PYTHON:-python3}
if command -v "$PY" >/dev/null 2>&1; then
    printf 'locus\tfamily\tmono\tend\tcons_start\tcons_end\tcovered\tfrac\tpartial\tpct_id\tstrand\tedge_gap\n' > "$W/terminal_monomers.tsv"
    printf 'FA\t146\nFB\t100\nFC\t50\n' > "$TMP/families.tsv"
    CONFIDENT_DELTA=10 "$PY" "$ROOT/satpipe/stages/05_subfamily_terminal_boundary/build_subfam_term_json.py" "$W" "$TMP/families.tsv" > "$TMP/json.log" 2>&1 \
        || { cat "$TMP/json.log"; echo "FAIL build_subfam_term_json.py exits 0"; exit 1; }
    J=$(cat "$W/subfamily.json")
    check "subfamily.json: L1 A2, confident (20 points)"   [ "${J#*'"L1":{"major":"FA","subfamily":"A2","pid":null,"bits":60.0,"second":"A1","delta":20.0,"confident":true}'}" != "$J" ]
    check "subfamily.json: L2 tie, not confident"          [ "${J#*'"L2":{"major":"FA","subfamily":"A3","pid":null,"bits":50.0,"second":"A1","delta":0.0,"confident":false}'}" != "$J" ]
    check "subfamily.json: L5 no call"                     [ "${J#*'"L5":{"major":"FA","subfamily":null,'}" != "$J" ]
else
    echo "skip     build_subfam_term_json.py ($PY not found)"
fi

# ---- the method switch (lib/satpipe/subx.sh)
source "$ROOT/lib/satpipe/common.sh"
source "$ROOT/lib/satpipe/flat.sh"
source "$ROOT/lib/satpipe/config.sh"
source "$ROOT/lib/satpipe/subx.sh"
require_gawk
CFG_FLAT=$TMP/cfg.flat
{ flat_emit subfamily map 1 file; flat_emit subfamily/method str bins file; } > "$CFG_FLAT"
check "method bins"                                [ "$(_p8_subfamily_method)" = bins ]
{ flat_emit subfamily map 1 file; flat_emit subfamily/method str boundary default; } > "$CFG_FLAT"
check "method boundary (schema default)"           [ "$(_p8_subfamily_method)" = boundary ]
{ flat_emit subfamily map 0 file; } > "$CFG_FLAT"
check "method absent: boundary"                    [ "$(_p8_subfamily_method)" = boundary ]
echo "all $N checks passed"
