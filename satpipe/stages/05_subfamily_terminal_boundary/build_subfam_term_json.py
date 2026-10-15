#!/usr/bin/env python3
"""Build subfamily.json, terminal_monomers.json and terminal_hist.json for the viewer.

Usage: build_subfam_term_json.py <work_dir> <families.tsv>
  work_dir      holds locus_subfamily.tsv (assign_subfamily.sh) and terminal_monomers.tsv
                (terminal_monomers.sh); receives the three JSON files
  families.tsv  family  mono  [consensus.fa]
Environment: CONFIDENT_DELTA  bits margin over the second-best sub-subfamily for a confident call [10]
"""
import collections
import json
import os
import sys

if len(sys.argv) != 3:
    sys.exit("Usage: build_subfam_term_json.py <work_dir> <families.tsv>")
W = sys.argv[1]
CONFIDENT_DELTA = float(os.environ.get("CONFIDENT_DELTA", "10"))
mono_of = {}
with open(sys.argv[2]) as f:
    for line in f:
        if not line.strip() or line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        mono_of[p[0]] = int(p[1])

# --- subfamily.json: locus -> {subfamily, pid, bits, second, delta, confident} ---
subfam = {}
with open(f"{W}/locus_subfamily.tsv") as f:
    next(f)
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) < 9: continue
        lid, major, sub, pid, aln, bits, second, sbits, delta = p
        subfam[lid] = {
            "major": major,
            "subfamily": None if sub == "NA" else sub,
            "pid": None if pid == "NA" else float(pid),
            "bits": None if bits == "NA" else float(bits),
            "second": None if second in ("NA","") else second,
            "delta": None if delta == "NA" else float(delta),
            # confident if margin over 2nd-best >= CONFIDENT_DELTA bits (or no competitor)
            "confident": (delta == "NA") or (delta != "NA" and float(delta) >= CONFIDENT_DELTA),
        }
with open(f"{W}/subfamily.json","w") as f:
    json.dump(subfam, f, separators=(",",":"))
print(f"subfamily.json: {len(subfam)} loci, {os.path.getsize(W+'/subfamily.json'):,} bytes")

# --- terminal_monomers.json: locus -> {start:{...}, end:{...}} ---
# fields per terminus: cons_start, cons_end, covered, frac, partial, pct_id, strand, edge_gap
term = {}
with open(f"{W}/terminal_monomers.tsv") as f:
    next(f)
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) < 10: continue
        lid, subt, mono, end, cs, ce, cov, frac, partial, pid = p[:10]
        strand = p[10] if len(p) > 10 else None      # written by terminal_monomers.sh since 2026-10-15
        edge_gap = (None if len(p) < 12 or p[11] == "NA" else int(p[11]))
        rec = {
            "strand": strand,
            "edge_gap": edge_gap,
            "cons_start": None if cs=="NA" else int(cs),
            "cons_end": None if ce=="NA" else int(ce),
            "covered": None if cov=="NA" else int(cov),
            "frac": None if frac=="NA" else float(frac),
            "partial": partial=="yes",
            "pct_id": None if pid=="NA" else float(pid),
            "mono": int(mono),
        }
        term.setdefault(lid, {})
        term[lid]["start" if end=="5p" else "end"] = rec
with open(f"{W}/terminal_monomers.json","w") as f:
    json.dump(term, f, separators=(",",":"))
print(f"terminal_monomers.json: {len(term)} loci, {os.path.getsize(W+'/terminal_monomers.json'):,} bytes")

# --- aggregate: covered-length histogram per family per end, for the viewer chart ---
agg = {}
for lid, ends in term.items():
    subt = subfam.get(lid,{}).get("major","?")
    for which in ("start","end"):
        r = ends.get(which)
        if not r or r["covered"] is None: continue
        agg.setdefault(subt, {}).setdefault(which, collections.Counter())
        # bin covered length into the monomer length range
        agg[subt][which][r["covered"]] += 1
# convert counters to arrays indexed by covered bp (0..mono)
out = {}
for subt, d in agg.items():
    if subt not in mono_of:
        print(f"  WARNING: family {subt!r} not in families.tsv; histogram skipped", file=sys.stderr)
        continue
    m = mono_of[subt]
    out[subt] = {"mono": m, "start":[0]*(m+2), "end":[0]*(m+2)}
    for which in ("start","end"):
        for cov,n in d.get(which, {}).items():
            if 0 <= cov <= m+1: out[subt][which][cov] += n
with open(f"{W}/terminal_hist.json","w") as f:
    json.dump(out, f, separators=(",",":"))
print(f"terminal_hist.json written")
for subt in sorted(out):
    s = out[subt]["start"]; e = out[subt]["end"]
    full_s = sum(s[int(out[subt]['mono']*0.9):])
    part_s = sum(s[:int(out[subt]['mono']*0.9)])
    print(f"  {subt}: 5' full={full_s} partial={part_s}")
