#!/usr/bin/env python3
"""
Comprehensive parental-subgenome analysis for discriminating CLsat sub-subfamilies.
For each discriminating subfamily, restricted to ANCHORABLE loci (>=1 unique flank):
  - locus counts per assembly
  - total monomer counts per assembly (n_monomers summed) -- caveated as collapse-sensitive
  - cross-assembly flank orthology: for each locus, which assemblies share an orthologous flank
"""
import json
from collections import defaultdict, Counter

FA="/staging/Sergei/Lacertidae/Darevskia/Satellites/manuscript/CLsat_full_clusters_6genome/flank_analysis"

# load locus master + subfamily + flank class
sp={}; subfam={}; nmono={}
for line in open(f"{FA}/locus_master.tsv"):
    p=line.rstrip("\n").split("\t")
    if len(p)<7: continue
    sp[p[0]]=p[4]; nmono[p[0]]=0
# n_monomers from flanks.json
flanks=json.load(open(f"{FA}/flanks.json"))
for loc,v in flanks.items():
    nmono[loc]=v.get('n_monomers',0)
# sub-subfamily
for line in open(f"{FA}/locus_subfamily.tsv"):
    p=line.rstrip("\n").split("\t")
    if p[0]=="locus" or len(p)<3: continue
    subfam[p[0]]=p[2]
# anchorable
anchor={}
for line in open(f"{FA}/locus_flank_class.tsv"):
    p=line.rstrip("\n").split("\t")
    if p[0]=="locus": continue
    anchor[p[0]]=(p[7]=="True")

# orthology
o=json.load(open(f"{FA}/flank_orthologs.json"))
THR="80"
# build per-locus set of assemblies it links to (union LF/RF), only counting unique-flank links
locus_links=defaultdict(set)
for side in ("LF","RF"):
    for q,hits in o[side][THR].items():
        for h in hits:
            if h["sp"]!=sp.get(q): locus_links[q].add(h["sp"])

SPECIES=["arm","mix","nai","unm","unp","val"]
PARENTAL=["nai","unm","unp","val"]

# discriminating CLsat sub-subfamilies
DISC=["CLsat1.1","CLsat1.2","CLsat3.1","CLsat3.2","CLsat1.5","CLsat4"]

print("="*70)
print("DISCRIMINATING CLsat SUB-SUBFAMILIES (anchorable loci only)")
print("="*70)
for d in DISC:
    loci=[l for l in subfam if subfam[l]==d and anchor.get(l,False)]
    all_loci=[l for l in subfam if subfam[l]==d]
    print(f"\n### {d}: {len(loci)} anchorable / {len(all_loci)} total loci")
    # per assembly: loci + monomers
    print(f"{'assembly':8} {'loci':>5} {'monomers':>9}  {'links_to (orthologous assemblies)'}")
    for s in SPECIES:
        sl=[l for l in loci if sp[l]==s]
        if not sl: continue
        mono=sum(nmono[l] for l in sl)
        # which assemblies do these loci link to
        link_counter=Counter()
        for l in sl:
            for tgt in locus_links.get(l,set()): link_counter[tgt]+=1
        links=", ".join(f"{a}:{n}" for a,n in sorted(link_counter.items(),key=lambda x:-x[1]))
        print(f"{s:8} {len(sl):>5} {mono:>9}  {links}")
