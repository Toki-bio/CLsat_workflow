#!/usr/bin/env bash
# diff_viewer.sh DIR_A DIR_B   compare two viewer directories (PORT_SPEC 4 viewer row, 7): same file set; every file
# byte-identical except *.gz, whose decompressed payloads must be identical (gzip member headers differ: Python's GzipFile
# stores the file name, `gzip -n` does not). Prints one line per difference and a summary; exit 1 on any difference.
set -uo pipefail
export LC_ALL=C
A=${1:?dir_a}; B=${2:?dir_b}
bad=0; n=0; nz=0
la=$(cd "$A" && find . -type f ! -name screenshot.png | LC_ALL=C sort); lb=$(cd "$B" && find . -type f ! -name screenshot.png | LC_ALL=C sort)
if [ "$la" != "$lb" ]; then echo "file sets differ:"; diff <(printf '%s\n' "$la") <(printf '%s\n' "$lb") | sed 's/^/  /'; bad=1; fi
while IFS= read -r f; do
    [ -n "$f" ] || continue
    [ -f "$B/$f" ] || continue
    n=$((n + 1))
    case $f in
        *.gz) nz=$((nz + 1))
              cmp -s <(gzip -dc < "$A/$f") <(gzip -dc < "$B/$f") || { echo "DIFF payload $f"; bad=1; } ;;
        *)    cmp -s "$A/$f" "$B/$f" || { echo "DIFF $f"; bad=1; } ;;
    esac
done <<< "$la"
echo "compared $n files ($nz gzip payloads): $([ $bad = 0 ] && echo IDENTICAL || echo DIFFERENT)"
exit $bad
