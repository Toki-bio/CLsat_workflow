#!/usr/bin/env python3
import json
from collections import defaultdict
W = "/staging/Sergei/Lacertidae/Darevskia/Satellites/manuscript/CLsat_full_clusters_6genome/flank_analysis"

# locus -> species
sp = {}
with open(f"{W}/locus_master.tsv") as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) < 5: continue
        sp[p[0]] = p[4]

o = json.load(open(f"{W}/flank_orthologs.json"))
SPECIES = ["arm","mix","nai","unm","unp","val"]
THR = "80"

# Count an orthologous LOCUS PAIR once if either flank links them (union of LF/RF).
# Build undirected pairs at species level.
pair_links = defaultdict(set)   # (spA,spB) -> set of frozenset(locusA,locusB)
for side in ("LF","RF"):
    d = o[side][THR]
    for qlocus, hits in d.items():
        qsp = sp.get(qlocus)
        if not qsp: continue
        for h in hits:
            tlocus = h["locus"]; tsp = h["sp"]
            if tsp == qsp: continue
            key = tuple(sorted([qsp, tsp]))
            pair_links[key].add(frozenset([qlocus, tlocus]))

print("=== Cross-assembly orthologous locus pairs (>=80% id, >=200bp; union of LF/RF), unique locus pairs ===")
# matrix
hdr = "      " + "".join(f"{s:>6}" for s in SPECIES)
print(hdr)
for a in SPECIES:
    row = f"{a:>6}"
    for b in SPECIES:
        if a == b:
            row += f"{'-':>6}"
        else:
            key = tuple(sorted([a,b]))
            row += f"{len(pair_links[key]):>6}"
    print(row)

print("\n=== Highlighted comparisons ===")
def npairs(a,b): return len(pair_links[tuple(sorted([a,b]))])
print(f"unm <-> unp (two haplotypes of same individual): {npairs('unm','unp')}")
print(f"unp <-> nai (paternal-hap vs nairensis):         {npairs('unp','nai')}")
print(f"unp <-> val (paternal-hap vs valentini):         {npairs('unp','val')}")
print(f"unm <-> nai (maternal-hap vs nairensis):         {npairs('unm','nai')}")
print(f"unm <-> val (maternal-hap vs valentini):         {npairs('unm','val')}")
print(f"nai <-> val (the two parents):                   {npairs('nai','val')}")
print(f"arm <-> mix (arm parentage incl. mixta):         {npairs('arm','mix')}")
print(f"arm <-> val (arm parentage incl. valentini):     {npairs('arm','val')}")
