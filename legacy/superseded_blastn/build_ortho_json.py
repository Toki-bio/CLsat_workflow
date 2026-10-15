#!/usr/bin/env python3
"""Build flank orthology JSON from blast results, at several identity thresholds.
Cross-species best reciprocal-ish hits per side (LF/RF)."""
import json, os
from collections import defaultdict

WORK = "/staging/Sergei/Lacertidae/Darevskia/Satellites/manuscript/CLsat_full_clusters_6genome/flank_analysis"

THRESHOLDS = [70, 80, 90]   # percent identity
MIN_ALN = 200               # min alignment length (bp)

def species_of(qid):
    # qid format "L123|sp"
    return qid.split("|", 1)[1] if "|" in qid else "?"
def locus_of(qid):
    return qid.split("|", 1)[0]

def load_blast(path):
    # qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    hits = []
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 12: continue
            q, s = p[0], p[1]
            if q == s: continue
            pid = float(p[2]); alen = int(p[3]); bits = float(p[11])
            hits.append((q, s, pid, alen, bits))
    return hits

def build(side):
    hits = load_blast(f"{WORK}/{side}.blast.tsv")
    # For each threshold, keep cross-species hits passing filters.
    # Store best hit per (query_locus -> target_species): the top-bitscore target locus.
    result = {}  # threshold -> { locus -> [ {locus,sp,pid,aln,bits}, ... ] }
    for thr in THRESHOLDS:
        best = defaultdict(dict)  # qlocus -> { tsp -> (bits, dict) }
        for q, s, pid, alen, bits in hits:
            if pid < thr or alen < MIN_ALN: continue
            qsp, ssp = species_of(q), species_of(s)
            if qsp == ssp: continue   # cross-species only
            ql, sl = locus_of(q), locus_of(s)
            cur = best[ql].get(ssp)
            if cur is None or bits > cur[0]:
                best[ql][ssp] = (bits, {"locus": sl, "sp": ssp,
                                        "pid": round(pid,1), "aln": alen,
                                        "bits": round(bits,1)})
        out = {}
        for ql, spmap in best.items():
            arr = sorted((v[1] for v in spmap.values()),
                         key=lambda d: -d["bits"])
            out[ql] = arr
        result[str(thr)] = out
    return result

ortho = {"LF": build("lf"), "RF": build("rf"), "min_aln": MIN_ALN,
         "thresholds": THRESHOLDS}

with open(f"{WORK}/flank_orthologs.json", "w") as f:
    json.dump(ortho, f, separators=(",", ":"))
print(f"wrote flank_orthologs.json ({os.path.getsize(WORK+'/flank_orthologs.json'):,} bytes)")

for side in ("LF", "RF"):
    for thr in THRESHOLDS:
        d = ortho[side][str(thr)]
        n_loci = len(d)
        n_pairs = sum(len(v) for v in d.values())
        print(f"  {side} @>={thr}%: {n_loci} loci have cross-species ortholog(s), {n_pairs} total links")
