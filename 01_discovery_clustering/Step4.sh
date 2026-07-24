#!/usr/bin/env bash
# run_trf_and_split_simple.sh — run TRF, compute repeat coverage, split FASTA
# Requires: gawk (for asorti)

set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <input.fasta>"
  exit 1
fi

fa="$1"
command -v trf >/dev/null 2>&1 || { echo "Error: trf not in PATH"; exit 1; }
command -v gawk >/dev/null 2>&1 || { echo "Error: gawk is required (asorti)"; exit 1; }

# Run TRF
trf "$fa" 2 5 7 80 10 50 2000 -h
dat="${fa}.2.5.7.80.10.50.2000.dat"
[[ -f "$dat" ]] || { echo "Error: dat file not found"; exit 1; }

# Parse TRF .dat -> seq coverage (merge overlaps correctly)
gawk '
function flush(    i, cov, cs, ce, si, ei, idx) {
  if (seq=="") return
  if (n==0) { print seq, 0; return }

  # sort indices by start[] values
  asorti(start, idx)

  cov=0
  cs=start[idx[1]]; ce=end[idx[1]]
  for (i=2;i<=n;i++) {
    si=start[idx[i]]; ei=end[idx[i]]
    if (si <= ce+1) { if (ei>ce) ce=ei }    # merge overlap/adjacent
    else { cov += ce - cs + 1; cs=si; ce=ei }
  }
  cov += ce - cs + 1
  print seq, cov
}

# New sequence header in .dat
/^Sequence:/ { flush(); seq=$2; n=0; delete start; delete end }
/^[0-9]/     { n++; start[n]=$1; end[n]=$2 }
END { flush() }
' "$dat" > "${fa}.coverage.txt"

# Split FASTA using coverage threshold 50 bp
gawk -v covfile="${fa}.coverage.txt" -v repfa="${fa}.repeats.fas" -v cleanfa="${fa}.clean.fas" '
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
    mode = (cov[hdr] > 50 ? "rep" : "clean")
  }
  if (mode=="rep")   print >> repfa
  if (mode=="clean") print >> cleanfa
}
' "$fa"

echo "Created:"
echo "  ${fa}.clean.fas  (coverage ≤ 50)"
echo "  ${fa}.repeats.fas (coverage > 50)"
