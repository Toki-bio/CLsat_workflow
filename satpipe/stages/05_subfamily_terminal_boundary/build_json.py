#!/usr/bin/env python3
"""Build viewer JSON files (flanks.json, boundary_phase.json, boundary_phase_hist.json) from the
flank-analysis work directory.

Usage: build_json.py <work_dir>
  work_dir holds all_regions.fa, locus_master.tsv, manifest.tsv, locus_boundary_phase.tsv
Note: this builder's flanks.json carries no flank `class`; the viewer's flanks.json comes from
rebuild_ortho_segmented.py (see PIPELINE.md). Kept for the boundary-phase products.
"""
import json
import os
import sys

if len(sys.argv) != 2:
    sys.exit("Usage: build_json.py <work_dir>")
WORK = sys.argv[1]

# --- Load region FASTA into dict: header -> seq ---
seqs = {}
hdr = None
buf = []
with open(f"{WORK}/all_regions.fa") as f:
    for line in f:
        line = line.rstrip("\n")
        if line.startswith(">"):
            if hdr is not None:
                seqs[hdr] = "".join(buf)
            hdr = line[1:]
            buf = []
        else:
            buf.append(line)
    if hdr is not None:
        seqs[hdr] = "".join(buf)
print(f"loaded {len(seqs)} region sequences")

# --- locus_master: Lid scf start0 end species subtype MONO ---
locus = {}
with open(f"{WORK}/locus_master.tsv") as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) < 7: continue
        lid, scf, s0, e, sp, subt, mono = p[:7]
        locus[lid] = dict(scf=scf, start=int(s0), end=int(e), sp=sp,
                          subtype=subt, mono=int(mono),
                          length=int(e)-int(s0))

# --- manifest: locus kind sp scf rs re trunc explen ---
flanks = {lid: {"sp": v["sp"], "scf": v["scf"], "subtype": v["subtype"],
                "start": v["start"], "end": v["end"], "length": v["length"],
                "mono": v["mono"],
                "n_monomers": round(v["length"]/v["mono"]) if v["mono"] else None,
                "LF": None, "RF": None}
          for lid, v in locus.items()}

with open(f"{WORK}/manifest.tsv") as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) < 8: continue
        lid, kind, sp, scf, rs, re_, trunc, explen = p
        if kind not in ("LF", "RF"): continue
        if rs == "0":
            continue  # missing flank (off scaffold edge)
        if lid not in flanks:
            continue
        header = f"{sp}|{scf}:{rs}-{re_}"
        seq = seqs.get(header, "")
        flanks[lid][kind] = {
            "seq": seq, "start": int(rs), "end": int(re_),
            "truncated": trunc == "1", "len": len(seq)
        }

# Mark missing flanks explicitly (start at scaffold edge)
n_lf = sum(1 for v in flanks.values() if v["LF"])
n_rf = sum(1 for v in flanks.values() if v["RF"])
print(f"loci with LF: {n_lf}, with RF: {n_rf}, total loci: {len(flanks)}")

with open(f"{WORK}/flanks.json", "w") as f:
    json.dump(flanks, f, separators=(",", ":"))
print(f"wrote flanks.json ({os.path.getsize(WORK+'/flanks.json'):,} bytes)")

# --- boundary phase JSON ---
bphase = {}
with open(f"{WORK}/locus_boundary_phase.tsv") as f:
    next(f)  # header
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) < 8: continue
        lid, kind, subt, cc, alen, pid, strand, qlen = p
        bphase.setdefault(lid, {})
        bphase[lid][kind] = {
            "cons_coord": None if cc == "NA" else int(cc),
            "aln_len": None if alen == "NA" else int(alen),
            "pct_id": None if pid == "NA" else float(pid),
            "strand": strand,
        }
with open(f"{WORK}/boundary_phase.json", "w") as f:
    json.dump(bphase, f, separators=(",", ":"))
print(f"wrote boundary_phase.json ({os.path.getsize(WORK+'/boundary_phase.json'):,} bytes)")

# --- aggregate stats: start/end consensus coord histograms per subtype ---
agg = {}
for lid, kinds in bphase.items():
    if lid not in locus:
        continue
    subt = locus[lid]["subtype"]
    agg.setdefault(subt, {"start": [0]* (locus[lid]["mono"]+1),
                          "end": [0]*(locus[lid]["mono"]+1),
                          "mono": locus[lid]["mono"]})
    lb = kinds.get("LB", {}).get("cons_coord")
    rb = kinds.get("RB", {}).get("cons_coord")
    if lb is not None and 1 <= lb <= locus[lid]["mono"]:
        agg[subt]["start"][lb] += 1
    if rb is not None and 1 <= rb <= locus[lid]["mono"]:
        agg[subt]["end"][rb] += 1
with open(f"{WORK}/boundary_phase_hist.json", "w") as f:
    json.dump(agg, f, separators=(",", ":"))
print(f"wrote boundary_phase_hist.json")
for subt in sorted(agg):
    s = agg[subt]["start"]; e = agg[subt]["end"]
    print(f"  {subt}: start peak at coord {s.index(max(s))} ({max(s)} loci), end peak at coord {e.index(max(e))} ({max(e)} loci)")
