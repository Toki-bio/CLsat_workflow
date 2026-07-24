#!/usr/bin/env python3
import json, os
W = "/staging/Sergei/Lacertidae/Darevskia/Satellites/manuscript/CLsat_full_clusters_6genome/flank_analysis"

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
            # confident if margin over 2nd-best >= 10 bits (or no competitor)
            "confident": (delta == "NA") or (delta != "NA" and float(delta) >= 10),
        }
with open(f"{W}/subfamily.json","w") as f:
    json.dump(subfam, f, separators=(",",":"))
print(f"subfamily.json: {len(subfam)} loci, {os.path.getsize(W+'/subfamily.json'):,} bytes")

# --- terminal_monomers.json: locus -> {start:{...}, end:{...}} ---
# fields per terminus: cons_start, cons_end, covered, frac, partial, pct_id
term = {}
with open(f"{W}/terminal_monomers.tsv") as f:
    next(f)
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) < 10: continue
        lid, subt, mono, end, cs, ce, cov, frac, partial, pid = p
        rec = {
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

# --- aggregate: covered-length histogram per subtype per end, for the viewer chart ---
agg = {}
import collections
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
mono_of = {"CLsat1":146,"CLsat2":145,"CLsat3":146,"CLsat4":116}
for subt, d in agg.items():
    m = mono_of.get(subt,146)
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
