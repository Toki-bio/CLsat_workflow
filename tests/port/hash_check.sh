#!/usr/bin/env bash
# PORT_SPEC 5(d): for every <run>/*/DONE written by the Python driver,
#   (1) re-serialise the embedded params (JSON -> TYPED flat -> json_dumps) and compare sha256[:16] with params_hash;
#   (2) re-emit the whole record with json_indent1 and compare it byte for byte with the file (json.dump indent=1 layout).
# Usage: tests/port/hash_check.sh RUN_DIR [RUN_DIR...]      exit 0 only when every DONE passes both checks
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
ROOT=$(cd -- "$HERE/../.." && pwd)
source "$ROOT/lib/satpipe/common.sh"
source "$ROOT/lib/satpipe/flat.sh"
source "$ROOT/lib/satpipe/done.sh"
require_gawk
[ $# -ge 1 ] || die "usage: hash_check.sh RUN_DIR [RUN_DIR...]" 2
TMP=$(mktemp -d "${TMPDIR:-/tmp}/satpipe_wp0.XXXXXX")
trap 'rm -rf -- "$TMP"' EXIT
total=0; okh=0; oki=0
for run in "$@"; do
    for d in "$run"/*/DONE; do
        [ -f "$d" ] || continue
        total=$((total + 1))
        stored=$(pyjson_gawk 'BEGIN { RS = "^$" } { json_parse($0, "D"); printf "%s", JV["D" SUBSEP "params_hash"] }' "$d")
        json_to_flat "$d" params > "$TMP/p.flat"
        got=$(params_hash "$TMP/p.flat")
        if [ "$got" = "$stored" ]; then okh=$((okh + 1)); else printf 'HASH MISMATCH %s stored=%s bash=%s\n' "$d" "$stored" "$got"; fi
        pyjson_gawk 'BEGIN { RS = "^$" } { json_parse($0, "D"); printf "%s", json_indent1("D", 0) }' "$d" > "$TMP/re.json"
        if cmp -s "$TMP/re.json" "$d"; then oki=$((oki + 1)); else printf 'INDENT1 DIFFERS %s\n' "$d"; fi
    done
done
printf 'params_hash reproduced: %d/%d    indent=1 record bytes reproduced: %d/%d\n' "$okh" "$total" "$oki" "$total"
[ "$total" -gt 0 ] && [ "$okh" = "$total" ] && [ "$oki" = "$total" ]
