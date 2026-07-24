#!/usr/bin/env python3
"""
Classify each CLsat locus flank as: truncated (scaffold edge), satellite-like, or unique.
A locus is 'anchorable' only if it has >=1 UNIQUE flank.
Loci with both flanks truncated-or-satellite are excluded from orthology-based parental analysis.
"""
import json, subprocess, os, tempfile

FA="/staging/Sergei/Lacertidae/Darevskia/Satellites/manuscript/CLsat_full_clusters_6genome/flank_analysis"
D6="/staging/Sergei/Lacertidae/Darevskia/Satellites/manuscript/CLsat_full_clusters_6genome"

flanks=json.load(open(f"{FA}/flanks.json"))

# satellite consensus bank: CLsat1-4 + DarSat all
cons_files=[f"{D6}/CLsat{i}.cons.fa" for i in (1,2,3,4)]
darsat=f"/staging/Sergei/Lacertidae/Darevskia/Satellites/CLUSTERS/Darevskia/DarSat_all_consensus.fa"
bank=tempfile.NamedTemporaryFile(mode='w',suffix='.fa',delete=False,dir=FA)
for cf in cons_files:
    if os.path.exists(cf):
        name=os.path.basename(cf).split('.')[0]
        seq="".join(l.strip() for l in open(cf) if not l.startswith(">"))
        bank.write(f">{name}\n{seq}{seq}\n")  # doubled
if os.path.exists(darsat):
    nm=None;sq=""
    for l in open(darsat):
        if l.startswith(">"):
            if nm: bank.write(f">{nm}\n{sq}{sq}\n")
            nm=l[1:].strip();sq=""
        else: sq+=l.strip()
    if nm: bank.write(f">{nm}\n{sq}{sq}\n")
bank.close()

# write all flank seqs to query fasta
qf=tempfile.NamedTemporaryFile(mode='w',suffix='.fa',delete=False,dir=FA)
flank_meta={}  # qid -> (locus, side, truncated, length)
for loc,v in flanks.items():
    for side in ('LF','RF'):
        fl=v[side]
        seq=fl.get('seq','') if isinstance(fl,dict) else ''
        trunc=fl.get('truncated',False) if isinstance(fl,dict) else False
        qid=f"{loc}__{side}"
        flank_meta[qid]=(loc,side,trunc,len(seq))
        if len(seq)>=50:
            qf.write(f">{qid}\n{seq}\n")
qf.close()

# ssearch36 flanks vs satellite bank; a flank is 'satellite-like' if it has a strong hit
# covering a large fraction of the flank
out=subprocess.run(["ssearch36","-n","-m","8C","-E","1e-3",qf.name,bank.name],
                   capture_output=True,text=True)
sat_hit={}  # qid -> best (bits, alnlen)
for line in out.stdout.splitlines():
    if line.startswith("#") or not line.strip(): continue
    p=line.split("\t")
    if len(p)<12: continue
    qid=p[0]; alen=int(p[3]); bits=float(p[11])
    cur=sat_hit.get(qid,(0,0))
    if bits>cur[0]: sat_hit[qid]=(bits,alen)

# classify
def classify(qid):
    loc,side,trunc,length=flank_meta[qid]
    if length<50 or trunc: return "truncated"
    bits,alen=sat_hit.get(qid,(0,0))
    # satellite-like if >=60% of the 500bp flank aligns to a satellite consensus at good score
    if alen>=300 and bits>=80: return "satellite"
    return "unique"

# per-locus result
rows=[]
for loc,v in flanks.items():
    lc=classify(f"{loc}__LF")
    rc=classify(f"{loc}__RF")
    n_unique=sum(1 for c in (lc,rc) if c=="unique")
    anchorable = n_unique>=1
    rows.append((loc,v['sp'],v['subtype'],v.get('n_monomers',0),lc,rc,n_unique,anchorable))

with open(f"{FA}/locus_flank_class.tsv","w") as o:
    o.write("locus\tsp\tsubtype\tn_monomers\tLF_class\tRF_class\tn_unique_flanks\tanchorable\n")
    for r in rows:
        o.write("\t".join(str(x) for x in r)+"\n")

# summary
from collections import Counter
total=len(rows)
anchorable=sum(1 for r in rows if r[7])
both_sat=sum(1 for r in rows if r[6]==0)
print(f"total loci: {total}")
print(f"anchorable (>=1 unique flank): {anchorable}")
print(f"EXCLUDED (both flanks truncated-or-satellite): {both_sat}")
print()
print("flank class combos (LF,RF):")
for combo,n in Counter((r[4],r[5]) for r in rows).most_common():
    print(f"  {combo}: {n}")
print()
print("anchorable loci per species:")
for sp,n in Counter(r[1] for r in rows if r[7]).most_common():
    print(f"  {sp}: {n}")
print()
print("EXCLUDED loci per species:")
for sp,n in Counter(r[1] for r in rows if not r[7]).most_common():
    print(f"  {sp}: {n}")

os.unlink(bank.name); os.unlink(qf.name)
