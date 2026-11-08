#!/usr/bin/env bash
# WP9: anchor phasing in bash/gawk (lib/satpipe/anchor.awk, anchor.sh).
#  A  the 21 cases of tests/test_anchor.py (6 tests, 16 of them parametrised rotations), through `satpipe anchor` (anchor_cli)
#  B  stage 21_anchor (st_21_anchor/plan_21_anchor with a minimal STUB of the WP3/WP5 driver functions, defined below)
#     against golden output of Python's st_anchor + plan (tests/sh/golden/anchor/, made at tag python-driver-final by
#     tests/port/anchor_stage_oracle.py): files, say line, manifest, outputs, plan hash, errors
#  C  the manuscript's 9 sub-subfamily consensuses re-rotate to fig1b_rotated.fasta (HANDOFF.md:216-219)
# Usage: tests/sh/test_anchor.sh      env: MANUSCRIPT_DIR (default C:/work/CLsat/manuscript)
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
ROOT=$(cd -- "$HERE/../.." && pwd)
source "$ROOT/lib/satpipe/common.sh"
source "$ROOT/lib/satpipe/flat.sh"
source "$ROOT/lib/satpipe/done.sh"
source "$ROOT/lib/satpipe/anchor.sh"
require_gawk
GOLDA=$ROOT/tests/sh/golden/anchor
MS=${MANUSCRIPT_DIR:-C:/work/CLsat/manuscript}
TMP=$(mktemp -d "${TMPDIR:-/tmp}/satpipe_wp9t.XXXXXX")
trap 'rm -rf -- "$TMP"' EXIT
CONS=$ROOT/satpipe/stages/02_seed_expansion_scan/consensus
if command -v cygpath >/dev/null 2>&1; then CONS=$(cygpath -m "$CONS"); TMP=$(cygpath -m "$TMP"); fi   # Windows python needs C:/ paths
MOTIF=AAAGAAG
N=0
check() { local label=$1; shift; N=$((N + 1)); if "$@"; then printf 'ok   %3d %s\n' "$N" "$label"; else printf 'FAIL %3d %s\n' "$N" "$label"; exit 1; fi; }
not() { ! "$@"; }
lib() { "$GAWK" -f "$PYJSON" -f "$SATPIPE_LIB/anchor.awk" -e "$@"; }
cons_seq() { P=$CONS/$1.cons.fa lib 'BEGIN { anc_read_fasta(ENVIRON["P"], R); printf "%s", anc_join_values(R) }'; }
seq_of() { P=$1 lib 'BEGIN { anc_read_fasta(ENVIRON["P"], R); printf "%s", anc_join_values(R) }'; }
quiet() { "$@" >/dev/null 2>&1; }

# config fixture: JSON (as the validated config data) -> TYPED flat (the input anchor_cli gets after run_check)
mkcfg() {   # mkcfg NAME MOTIF_JSON [EXTRA_JSON]
    local f fams="" sep=""
    for f in CLsat1 CLsat2 CLsat3 CLsat4; do fams+="$sep\"$f\": {\"consensus\": \"$CONS/$f.cons.fa\", \"colour\": null}"; sep=", "; done
    printf '{"run_name": "t", "families": {%s}, "anchor": {"family": "CLsat1", "motif": %s, "independent": false, "max_motif_mismatch": 1, "min_identity": 0.6}%s}\n' \
        "$fams" "$2" "${3:-}" > "$TMP/$1.json"
    json_to_flat "$TMP/$1.json" > "$TMP/$1.flat"
}
mkcfg base '"AAAGAAG"'
mkcfg nomotif '""'
mkcfg absent '"CCCCCCCCCC"'
FLAT=$TMP/base.flat
fa() { printf '>%s\n%s\n' "$2" "$3" > "$1"; }

echo "== A: tests/test_anchor.py"
FRAME=$(S=$(cons_seq CLsat1) lib 'BEGIN { printf "%s", anc_frame(ENVIRON["S"], "AAAGAAG", 1) }')
in_frame_all() {
    [ "${FRAME:0:7}" = "$MOTIF" ] || return 1
    local f; for f in CLsat1 CLsat2 CLsat3 CLsat4; do quiet anchor_cli "$FLAT" --check "$CONS/$f.cons.fa" || return 1; done
}
check "test_shipped_consensuses_are_in_frame" in_frame_all
rot_case() {   # test_any_rotation_is_restored[fam-shift]
    local fam=$1 shift=$2 s k r d=$TMP/rot_$1_$2
    s=$(cons_seq "$fam"); k=$((shift % ${#s})); r=${s:k}${s:0:k}
    fa "$d.fa" x "$r"
    not quiet anchor_cli "$FLAT" --check "$d.fa" || return 1                 # out of frame is detected
    quiet anchor_cli "$FLAT" "$d.fa" -o "$d.out.fa" || return 1
    [ "$(seq_of "$d.out.fa")" = "$s" ] || return 1
    quiet anchor_cli "$FLAT" --check "$d.out.fa"
}
for fam in CLsat1 CLsat2 CLsat3 CLsat4; do for sh in 1 36 70 111; do check "test_any_rotation_is_restored[$fam-$sh]" rot_case "$fam" "$sh"; done; done
gapped() { local s; s=$(cons_seq CLsat2); fa "$TMP/g.fa" g "${s:0:10}--${s:10}"; quiet anchor_cli "$FLAT" --check "$TMP/g.fa"; }
check "test_gaps_ignored_in_check" gapped
JUNK=CAGATTTTCATATTATGCAGAAAATCTACTTCGCCTGATACGAGTCGGTTATCTTCGGATACTGTATAGTCCCACCTGGTGATCCTATGCTTGTGAGTACCCAGAAAATAGCGACGGACCGCGGTGTTAAGTGTCGAGCTACATCACTTC   # random.seed(1), 150 x choice('ACGT')
junk() { fa "$TMP/junk.fa" junk "$JUNK"; local e; e=$(anchor_cli "$FLAT" "$TMP/junk.fa" 2>&1 >/dev/null) && return 1; [[ $e == "satpipe anchor: anchor phasing failed:"* ]]; }
check "test_unplaceable_sequence_fails" junk
absent() { local e; e=$(anchor_cli "$TMP/absent.flat" --check "$CONS/CLsat1.cons.fa" 2>&1) && return 1
           [ "$e" = "satpipe anchor: anchor motif CCCCCCCCCC not found in the anchor consensus (best 6 mismatches)" ] || { echo "$e"; return 1; }; }
check "test_motif_absent_from_anchor_is_an_error" absent
nomotif() {
    local a s f; a=$(cons_seq CLsat1); s=$(cons_seq CLsat3)
    f=$(S=$a lib 'BEGIN { printf "%s", anc_frame(ENVIRON["S"], "", 1) }'); [ "$f" = "$a" ] || return 1
    fa "$TMP/nm.fa" x "${s:50}${s:0:50}"
    quiet anchor_cli "$TMP/nomotif.flat" "$TMP/nm.fa" -o "$TMP/nm.out.fa" && quiet anchor_cli "$TMP/nomotif.flat" --check "$TMP/nm.out.fa"
}
check "test_no_motif_uses_anchor_as_given" nomotif
echo "A: $N/21 pytest cases pass"

echo "== B: stage 21_anchor vs Python st_anchor (driver functions STUBBED here: cfg, cfg_type, cfg_keys, cfg_subtree, emit_inputs, emit_outputs)"
cfg() { flat_get "$CFG_FLAT" "$1"; }
cfg_type() { flat_type "$CFG_FLAT" "$1"; }
cfg_keys() { flat_keys "$CFG_FLAT" "$1"; }
cfg_subtree() { flat_subtree "$CFG_FLAT" "$1" "$2" >> "$PLAN_PARAMS"; }
emit_inputs() { printf '%s\n' "$@" >> "$PLAN_INPUTS"; }
emit_outputs() { printf '%s\n' "$@" >> "$STAGE_OUTPUTS"; }
manifest_set() { printf '%s' "$2" > "$RUN/manifest_anchor.json"; }

# sub-subfamily inputs: the manuscript FASTA when present, otherwise the four family consensuses renamed
SUBS=$MS/subsubfamily_consensuses.fasta
mkdir -p "$TMP/subdir"
if [ -f "$SUBS" ]; then tr -d '\r' < "$SUBS" > "$TMP/subs.fa"; else
    for f in CLsat1 CLsat2 CLsat3 CLsat4; do printf '>%s.1\n%s\n' "$f" "$(cons_seq $f)"; done > "$TMP/subs.fa"; fi
NAMES=$(grep '^>' "$TMP/subs.fa" | cut -c2- | tr '\n' ' ')
members_json() {   # members grouped by the family prefix of each record name, config order
    local n fam out="" cur="" sep=""
    for n in $NAMES; do fam=${n%%.*}
        if [ "$fam" != "$cur" ]; then [ -n "$cur" ] && out+="]"; out+="$sep\"$fam\": [\"$n\""; cur=$fam; sep=", "; else out+=", \"$n\""; fi
    done
    printf '{%s]}' "$out"
}
MEM=$(members_json)
# directory input: two files + an .anchored. copy and a non-FASTA file that fasta_files() must skip, + a dotfile it keeps
gawk 'NR % 4 == 1 || NR % 4 == 2' "$TMP/subs.fa" > "$TMP/subdir/a.fa"; gawk 'NR % 4 == 3 || NR % 4 == 0' "$TMP/subs.fa" > "$TMP/subdir/B.FASTA"
cp "$TMP/subs.fa" "$TMP/subdir/x.anchored.fa"; echo junk > "$TMP/subdir/notes.txt"; : > "$TMP/subdir/.empty.fa"
fa "$TMP/dup.fa" CLsat1.1 ACGT

stage_case_json() {   # like stage_case on an existing $TMP/st_NAME.json
    local name=$1 c=$TMP/st_$1
    json_to_flat "$TMP/st_$name.json" > "$TMP/st_$name.flat"; CFG_FLAT=$TMP/st_$name.flat; mkdir -p "$c/sh"
    (   RUN=$c/sh PLAN_INPUTS=$c/inputs PLAN_PARAMS=$c/params STAGE_OUTPUTS=$c/outputs; : > "$PLAN_INPUTS"; : > "$PLAN_PARAMS"; : > "$STAGE_OUTPUTS"
        plan_21_anchor
        printf 'PLAN_HASH\t%s\n' "$(params_hash "$PLAN_PARAMS")"
        sed 's/^/INPUT\t/' "$PLAN_INPUTS"
        msg=$( { st_21_anchor > "$c/say.txt"; } 2>&1 ) || { M=${msg#satpipe run: } "$GAWK" -f "$PYJSON" -e 'BEGIN { printf "ERROR\t%s\n", flat_esc(ENVIRON["M"]) }'; exit 0; }
        sed 's/^/SAY\t/' "$c/say.txt"
        printf 'MANIFEST\t'; json_to_flat "$RUN/manifest_anchor.json" > "$c/m.flat"; flat_json "$c/m.flat"; echo
        sed 's/^/OUTPUT\t/' "$STAGE_OUTPUTS"
    ) | sed "s#$c/sh#@RUN@#g; s#$TMP#@TMP@#g; s#$CONS#@CONS@#g; s#^PLAN_HASH	.*#PLAN_HASH	@HASH@#" > "$c/sh.out"   # params hold machine paths
    diff "$GOLDA/$name.out" "$c/sh.out" || return 1
    if [ -d "$GOLDA/$name.anchored" ]; then diff -r --strip-trailing-cr "$GOLDA/$name.anchored" "$c/sh/anchored" || return 1; fi
    return 0
}
stage_case() { mkcfg "st_$1" "$2" "$3"; stage_case_json "$1"; }   # stage_case NAME MOTIF_JSON EXTRA_JSON
check "stage: Darevskia frame, 4 families + sub-subfamilies (one file)" stage_case dar '"AAAGAAG"' ", \"subfamily\": {\"members\": $MEM, \"consensuses\": \"$TMP/subs.fa\"}"
check "stage: sub-subfamilies from a directory (.anchored./non-FASTA skipped, dotfile kept)" stage_case dir '"AAAGAAG"' ", \"subfamily\": {\"members\": $MEM, \"consensuses\": \"$TMP/subdir\"}"
check "stage: independent families (own frame, subs onto their family)" stage_case indep '""' ", \"subfamily\": {\"members\": $MEM, \"consensuses\": \"$TMP/subs.fa\"}"
indep_true() { sed 's/"independent": false/"independent": true/' "$TMP/st_indep.json" > "$TMP/st_indep2.json"; stage_case_json indep2; }
check "stage: independent: true" indep_true
errcase() { stage_case "$@" && grep -q '^ERROR' "$GOLDA/$1.out"; }
check "stage error: subfamily.consensuses not found" errcase miss '"AAAGAAG"' ", \"subfamily\": {\"members\": null, \"consensuses\": \"$TMP/nowhere.fa\"}"
cat "$TMP/subs.fa" "$TMP/dup.fa" > "$TMP/subdir2.fa"; mkdir -p "$TMP/dupdir"; cp "$TMP/subs.fa" "$TMP/dupdir/1.fa"; cp "$TMP/dup.fa" "$TMP/dupdir/2.fa"
check "stage error: record appears twice (one file: ValueError text)" errcase dup1 '"AAAGAAG"' ", \"subfamily\": {\"members\": null, \"consensuses\": \"$TMP/subdir2.fa\"}"
check "stage error: record appears twice across files (StageError)" errcase dup2 '"AAAGAAG"' ", \"subfamily\": {\"members\": null, \"consensuses\": \"$TMP/dupdir\"}"
check "stage error: listed sub-subfamily without a record" errcase lost '"AAAGAAG"' ", \"subfamily\": {\"members\": {\"CLsat1\": [\"CLsat1.1\", \"CLsat9.9\", \"X\"]}, \"consensuses\": \"$TMP/subs.fa\"}"
check "stage error: motif absent from the anchor (AnchorError)" errcase nomo '"CCCCCCCCCC"' ""
check "stage error: ambiguous motif (AnchorError)" errcase amb '"TTTG"' ""
fa "$TMP/junkfam.fa" junk "$JUNK"
junkfam() { mkcfg st_jf '"AAAGAAG"' ""; sed -i "s#$CONS/CLsat3.cons.fa#$TMP/junkfam.fa#" "$TMP/st_jf.json"; stage_case_json jf && grep -q '^ERROR' "$GOLDA/jf.out"; }
check "stage error: a family cannot be placed (AnchorError lists it)" junkfam
nofam() { mkcfg st_nf '""' ""; sed -i 's/"family": "CLsat1"/"family": null/' "$TMP/st_nf.json"; stage_case_json nf; }
check "stage: anchor.family null -> first family, no motif" nofam
noblock() { printf '{"run_name": "t", "families": {"CLsat2": {"consensus": "%s"}, "CLsat1": {"consensus": "%s"}}}\n' "$CONS/CLsat2.cons.fa" "$CONS/CLsat1.cons.fa" > "$TMP/st_nb.json"; stage_case_json nb; }
check "stage: no anchor/subfamily blocks at all (defaults from .get)" noblock

echo "== C: manuscript sub-subfamily consensuses (HANDOFF.md:216-219)"
if [ -f "$MS/subsubfamily_consensuses.fasta" ] && [ -f "$MS/fig1b_rotated.fasta" ]; then
    ms() {
        tr -d '\r' < "$MS/subsubfamily_consensuses.fasta" > "$TMP/ms_in.fa"; tr -d '\r' < "$MS/fig1b_rotated.fasta" > "$TMP/ms_ref.fa"
        anchor_cli "$FLAT" "$TMP/ms_in.fa" -o "$TMP/ms_sh.fa" >/dev/null || return 1
        # the reference may be line-wrapped: compare record by record on unwrapped sequences
        "$GAWK" '/^>/ { if (n) print n "\t" s; n = substr($1, 2); s = ""; next } { s = s $0 } END { print n "\t" s }' "$TMP/ms_ref.fa" > "$TMP/ms_ref.tsv"
        "$GAWK" '/^>/ { if (n) print n "\t" s; n = substr($1, 2); s = ""; next } { s = s $0 } END { print n "\t" s }' "$TMP/ms_sh.fa" > "$TMP/ms_sh.tsv"
        diff "$TMP/ms_ref.tsv" "$TMP/ms_sh.tsv" || return 1
        [ "$(grep -c '' "$TMP/ms_sh.tsv")" = 9 ] || return 1
        quiet anchor_cli "$FLAT" --check "$TMP/ms_ref.fa"
    }
    check "9 sub-subfamily consensuses re-rotate to exactly fig1b_rotated.fasta; fig1b passes --check" ms
else
    echo "SKIP C: subsubfamily_consensuses.fasta / fig1b_rotated.fasta not found under MANUSCRIPT_DIR=$MS"
fi
echo "all $N checks pass"
