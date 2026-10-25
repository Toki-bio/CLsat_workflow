#!/usr/bin/env bash
# Step 4 - run TRF on the motifs themselves, compute each motif's internal simple-repeat
# coverage, and split the FASTA into clean / repeat-dominated motifs.
# Requires: trf, gawk (PROCINFO["sorted_in"])
#
# Usage: Step4.sh <input.fasta>
# Environment (defaults = published Darevskia run):
#   TRF_PARAMS          TRF parameters "match mismatch delta PM PI minscore maxperiod"  [2 5 7 80 10 50 2000]
#   SIMPLE_REPEAT_COV   motifs whose internal TRF coverage exceeds this many bp are
#                       written to .repeats.fas instead of .clean.fas                 [50]

set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <input.fasta>"
  exit 1
fi

fa="$1"
TRF_PARAMS="${TRF_PARAMS:-2 5 7 80 10 50 2000}"
COV_MAX="${SIMPLE_REPEAT_COV:-50}"
command -v trf >/dev/null 2>&1 || { echo "Error: trf not in PATH"; exit 1; }
command -v gawk >/dev/null 2>&1 || { echo "Error: gawk is required (PROCINFO sorted_in)"; exit 1; }
[[ "$TRF_PARAMS" =~ ^[0-9]+(\ [0-9]+){6}$ ]] || { echo "Error: TRF_PARAMS must be 7 integers, got '$TRF_PARAMS'"; exit 1; }

# Run TRF. TRF's exit status is not 0 on success in every build (some return the number of
# sequences processed), so do not let `set -e` abort here; the .dat check below is the real test.
dat="${fa}.${TRF_PARAMS// /.}.dat"
rm -f "$dat"   # a .dat left by an earlier run would pass the existence check below (GLM audit #51)
# shellcheck disable=SC2086
rc=0; trf "$fa" $TRF_PARAMS -h || rc=$?
# exit >= 125 = killed / crashed: the .dat is partial (lower non-zero codes are normal for some TRF builds)
[[ "$rc" -lt 125 ]] || { echo "Error: TRF crashed on $fa (exit $rc)" >&2; exit 1; }
[[ -f "$dat" ]] || { echo "Error: dat file not found: $dat"; exit 1; }

# Parse TRF .dat -> seq coverage (merge overlaps correctly)
gawk '
function flush(    k, cov, cs, ce, si, ei, first) {
  if (seq=="") return
  if (n==0) { print seq, 0; return }

  # Visit intervals in ascending order of their START VALUE.
  # (The earlier asorti(start, idx) sorted the index strings "1","10","11",...,"2", not the
  #  starts, so for >= 10 TRF hits the merge below double-counted overlaps; fixed 2026-10-15.)
  PROCINFO["sorted_in"] = "@val_num_asc"
  cov=0; first=1
  for (k in start) {
    si=start[k]; ei=end[k]
    if (first) { cs=si; ce=ei; first=0; continue }
    if (si <= ce+1) { if (ei>ce) ce=ei }    # merge overlap/adjacent
    else { cov += ce - cs + 1; cs=si; ce=ei }
  }
  cov += ce - cs + 1
  PROCINFO["sorted_in"] = "@unsorted"
  print seq, cov
}

# New sequence header in .dat
/^Sequence:/ { flush(); seq=$2; n=0; delete start; delete end }
/^[0-9]/     { n++; start[n]=$1; end[n]=$2 }
END { flush() }
' "$dat" > "${fa}.coverage.txt"

# Split FASTA using the coverage threshold
# Fresh outputs on every run (the awk below appends).
rm -f "${fa}.repeats.fas" "${fa}.clean.fas"
gawk -v covfile="${fa}.coverage.txt" -v repfa="${fa}.repeats.fas" -v cleanfa="${fa}.clean.fas" -v covmax="$COV_MAX" '
BEGIN {
  FS="[ \t]+"
  while ((getline line < covfile) > 0) {
    split(line,a,FS); cov[a[1]] = a[2]+0
  }
  close(covfile)
}
{
  if ($0 ~ /^>/) {
    hdr = substr($0,2)
    mode = (cov[hdr] > covmax ? "rep" : "clean")
  }
  if (mode=="rep")   print >> repfa
  if (mode=="clean") print >> cleanfa
}
' "$fa"
# both outputs always exist: a genome whose motifs are all internal repeats has an EMPTY clean set, it must not
# vanish from the later stages without a word (GLM audit 2026-10-24)
touch "${fa}.repeats.fas" "${fa}.clean.fas"
[ -s "${fa}.clean.fas" ] || echo "WARNING: ${fa}: no motif passed the internal-repeat filter (clean set is empty)" >&2

echo "Created:"
echo "  ${fa}.clean.fas  (coverage <= ${COV_MAX})"
echo "  ${fa}.repeats.fas (coverage > ${COV_MAX})"
