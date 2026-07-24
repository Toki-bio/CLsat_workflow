#!/usr/bin/env bash
set -euo pipefail

# arbitration_patched.sh -- CLsat locus family assignment pipeline
#
# Patches vs original:
#  D3 : pos-1 for 0-based BED start; CLsat1 window = 292 bp (MONO=146)
#  D3b: sorted() tiebreak in Python segmentation (deterministic)
#  D6 : three explicit branches: UNIQUE_SEED / SEED_CONFLICT+sub / NO_SEED
#       sub-loci and SEED_CONFLICT get delta scrutiny via DELTA_AMBIG
#       NO_SEED uses bitscore winner without phantom 0-delta comparison
#
# Output columns: locus_id family assigned_bits best_other best_other_bits
#                 delta flag provenance(col8)

BASE="/staging/tmp/clsat_v6/table1_540"
OUT="/staging/tmp/clsat_v6/table1_540"

BITS1="/staging/tmp/clsat_v6/pass3/scan_CLsat1.bits.tsv"
BITS2="/staging/tmp/clsat_v6/pass3/scan_CLsat2.bits.tsv"
BITS3="/staging/tmp/clsat_v6/pass3/scan_CLsat3.bits.tsv"
BITS4="/staging/tmp/clsat_v6/pass3/scan_CLsat4.bits.tsv"

SEED1="/staging/Sergei/Lacertidae/Darevskia/Satellites/motifs_output/clean/cross_clusters/selected/subclusters/seeds/CLsat1.bed"
SEED2="/staging/Sergei/Lacertidae/Darevskia/Satellites/motifs_output/clean/cross_clusters/selected/subclusters/seeds/CLsat2.bed"
SEED3="/staging/Sergei/Lacertidae/Darevskia/Satellites/motifs_output/clean/cross_clusters/selected/subclusters/seeds/CLsat3.bed"
SEED4="/staging/Sergei/Lacertidae/Darevskia/Satellites/motifs_output/clean/cross_clusters/selected/subclusters/seeds/CLsat4.bed"

if [[ -s "$OUT/MASTER.merged.bed" ]]; then
  MASTER_BED="$OUT/MASTER.merged.bed"; HAS_PROV=1
else
  MASTER_BED="$OUT/MASTER.merged.bed"; HAS_PROV=0
fi

DELTA=20        # UNIQUE_SEED contradiction threshold
DELTA_AMBIG=15  # CONFLICT/sub-locus/NO_SEED ambiguity threshold
SEG_BIN=1000; SEG_MIN_RUN=10; SEG_MIN_LOCUS=2000

die(){ echo "FATAL: $*" >&2; exit 1; }
check_file(){ [[ -s "$1" ]] || die "Missing: $1"; }
log(){ echo "[$(date +"%F %T")] $*"; }

log "Pre-flight (DELTA=${DELTA} DELTA_AMBIG=${DELTA_AMBIG} prov=${HAS_PROV})"
for f in "$BITS1" "$BITS2" "$BITS3" "$BITS4" "$SEED1" "$SEED2" "$SEED3" "$SEED4" "$MASTER_BED"; do
  check_file "$f"; done
command -v bedtools >/dev/null || die "bedtools not in PATH"
command -v python3  >/dev/null || die "python3 not in PATH"

##############################################################################
# D1
##############################################################################
log "D1: assigning locus IDs"

awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,"L"NR}' "$MASTER_BED" > "$OUT/MASTER.loci.bed"

if [[ "$HAS_PROV" == "1" ]]; then
  awk 'BEGIN{FS=OFS="\t"}{print "L"NR, ($5?$5:".")}' "$MASTER_BED" \
    > "$OUT/MASTER.locus_provenance.tsv"
else
  awk 'BEGIN{OFS="\t"}{print "L"NR, "."}' "$MASTER_BED" \
    > "$OUT/MASTER.locus_provenance.tsv"
fi
check_file "$OUT/MASTER.loci.bed"

##############################################################################
# D2
##############################################################################
log "D2: tagging seeds"

cat \
  <(awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,"CLsat1"}' "$SEED1") \
  <(awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,"CLsat2"}' "$SEED2") \
  <(awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,"CLsat3"}' "$SEED3") \
  <(awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,"CLsat4"}' "$SEED4") \
| sort -k1,1 -k2,2n > "$OUT/SEEDS.tagged.bed"
check_file "$OUT/SEEDS.tagged.bed"

##############################################################################
# D3 -- PATCHED: 0-based BED start; CLsat1=292bp
##############################################################################
log "D3: bits TSV -> BED (pos-1; CLsat1=292bp)"

awk 'BEGIN{FS=OFS="\t"} NR>1{print $1, $2-1, $2-1+292, $3}' "$BITS1" > "$OUT/CLsat1.bits.bed"
awk 'BEGIN{FS=OFS="\t"} NR>1{print $1, $2-1, $2-1+290, $3}' "$BITS2" > "$OUT/CLsat2.bits.bed"
awk 'BEGIN{FS=OFS="\t"} NR>1{print $1, $2-1, $2-1+292, $3}' "$BITS3" > "$OUT/CLsat3.bits.bed"
awk 'BEGIN{FS=OFS="\t"} NR>1{print $1, $2-1, $2-1+232, $3}' "$BITS4" > "$OUT/CLsat4.bits.bed"
for f in "$OUT"/CLsat*.bits.bed; do check_file "$f"; done

##############################################################################
# D3b -- PATCHED: deterministic sorted() tiebreak
##############################################################################
log "D3b: segmenting composite loci"

python3 - "$OUT/MASTER.loci.bed" "$BITS1" "$BITS2" "$BITS3" "$BITS4" \
  "$SEG_BIN" "$SEG_MIN_RUN" "$SEG_MIN_LOCUS" "$OUT/MASTER.loci.segmented.bed" <<'PYEOF'
import sys
from collections import defaultdict
loci_file=sys.argv[1]; bits_files=sys.argv[2:6]
families=["CLsat1","CLsat2","CLsat3","CLsat4"]
BIN=int(sys.argv[6]); MIN_RUN=int(sys.argv[7])
MIN_LOCUS=int(sys.argv[8]); out_file=sys.argv[9]
loci=[]
with open(loci_file) as f:
    for line in f:
        p=line.rstrip("\n").split("\t"); loci.append((p[0],int(p[1]),int(p[2]),p[3]))
print(f"  Loaded {len(loci)} loci",file=sys.stderr)
bins=defaultdict(lambda:defaultdict(lambda:defaultdict(float)))
for fi,bf in enumerate(bits_files):
    fam=families[fi]
    with open(bf) as f:
        f.readline()
        for line in f:
            p=line.rstrip("\n").split("\t"); scf,pos,bits=p[0],int(p[1]),float(p[2])
            b=pos//BIN
            if bits>bins[scf][b][fam]: bins[scf][b][fam]=bits
    print(f"  Loaded {bf} ({fam})",file=sys.stderr)
total_seg=0; n_split=0; out=open(out_file,"w")
for scf,start,end,lid in loci:
    bin_start=start//BIN; bin_end=(end-1)//BIN
    dom=[]
    for b in range(bin_start,bin_end+1):
        fb=bins[scf][b]
        # PATCHED: sorted() for deterministic alphabetical tiebreak
        dom.append(max(sorted(fb.keys()),key=fb.get) if fb else None)
    if not dom:
        out.write(f"{scf}\t{start}\t{end}\t{lid}\n"); total_seg+=1; continue
    segments=[]; seg_start_idx=0; cur_fam=dom[0]; i=0
    while i<len(dom):
        if dom[i] is not None and dom[i]!=cur_fam:
            new_fam=dom[i]; run=0; j=i
            while j<len(dom) and dom[j]==new_fam: run+=1; j+=1
            if run>=MIN_RUN:
                segments.append((seg_start_idx,i-1,cur_fam)); seg_start_idx=i; cur_fam=new_fam; i=j
            else: i+=1
        else: i+=1
    segments.append((seg_start_idx,len(dom)-1,cur_fam))
    if len(segments)==1:
        out.write(f"{scf}\t{start}\t{end}\t{lid}\n"); total_seg+=1
    else:
        n_split+=1
        for si,(sb,eb,fam) in enumerate(segments):
            s_bp=max(start,(bin_start+sb)*BIN); e_bp=min(end,(bin_start+eb+1)*BIN)
            if e_bp-s_bp<MIN_LOCUS: continue
            out.write(f"{scf}\t{s_bp}\t{e_bp}\t{lid}.{si+1}\n"); total_seg+=1
out.close()
print(f"  Original: {len(loci)}, split: {n_split}, output: {total_seg}",file=sys.stderr)
PYEOF

check_file "$OUT/MASTER.loci.segmented.bed"
LOCI_BED="$OUT/MASTER.loci.segmented.bed"

awk 'BEGIN{FS=OFS="\t"}
FNR==NR{ prov[$1]=$2; next }
{ lid=$4; parent=lid; sub(/\.[0-9]+$/,"",parent)
  p=(parent in prov?prov[parent]:(lid in prov?prov[lid]:"."))
  print lid,p }' \
"$OUT/MASTER.locus_provenance.tsv" "$LOCI_BED" \
> "$OUT/MASTER.locus_provenance.segmented.tsv"

##############################################################################
# D2b
##############################################################################
log "D2b: seeds -> segmented loci"

bedtools intersect -a "$LOCI_BED" -b "$OUT/SEEDS.tagged.bed" -wa -wb \
| awk 'BEGIN{FS=OFS="\t"}
  { key=$4 SUBSEP $8; seen[key]=1 }
  END{ for(k in seen){ split(k,a,SUBSEP); print a[1],a[2] } }' \
> "$OUT/MASTER.locus_seed_hits.tsv"
check_file "$OUT/MASTER.locus_seed_hits.tsv"

awk 'BEGIN{FS=OFS="\t"}
{ hits[$1]=(hits[$1]?hits[$1]","$2:$2); n[$1]++ }
END{
  print "locus_id","seed_call","seed_hits","seed_n"
  for(l in n){
    if(n[l]==1) print l,hits[l],hits[l],1
    else        print l,"SEED_CONFLICT",hits[l],n[l]
  }
}' "$OUT/MASTER.locus_seed_hits.tsv" > "$OUT/MASTER.locus_seeds.tsv"
check_file "$OUT/MASTER.locus_seeds.tsv"
log "  Seed calls: $(wc -l < "$OUT/MASTER.locus_seeds.tsv")"
log "  SEED_CONFLICT: $(grep -c SEED_CONFLICT "$OUT/MASTER.locus_seeds.tsv" || echo 0)"

##############################################################################
# D4
##############################################################################
log "D4: max bits per locus"

for i in 1 2 3 4; do
  bedtools intersect -a "$LOCI_BED" -b "$OUT/CLsat${i}.bits.bed" -wa -wb \
  | awk 'BEGIN{FS=OFS="\t"}
    { l=$4; b=$8; if(b>max[l]) max[l]=b }
    END{ for(l in max) print l,max[l] }' \
  > "$OUT/MASTER.CLsat${i}.maxbits.tsv"
  check_file "$OUT/MASTER.CLsat${i}.maxbits.tsv"
done

##############################################################################
# D5
##############################################################################
log "D5: locus bits matrix"

cut -f4 "$LOCI_BED" | awk '{print $1,0,0,0,0}' OFS="\t" > "$OUT/MASTER.locus_bits.tsv"

for i in 1 2 3 4; do
  col=$((i+1))
  awk -v C="$col" 'BEGIN{FS=OFS="\t"}
  FNR==NR{b[$1]=$2; next}
  { if($1 in b) $C=b[$1]; print }
  ' "$OUT/MASTER.CLsat${i}.maxbits.tsv" "$OUT/MASTER.locus_bits.tsv" \
  > "$OUT/tmp" && mv "$OUT/tmp" "$OUT/MASTER.locus_bits.tsv"
done

##############################################################################
# D6 -- PATCHED: three explicit branches
##############################################################################
log "D6: arbitration (three-branch)"

awk -v D="$DELTA" -v DA="$DELTA_AMBIG" 'BEGIN{FS=OFS="\t"}

FNR==NR{ if(FNR==1) next; seed[$1]=$2; next }

{
  l=$1
  bits["CLsat1"]=$2+0; bits["CLsat2"]=$3+0
  bits["CLsat3"]=$4+0; bits["CLsat4"]=$5+0
  s=seed[l]

  # Branch 1: SEED_CONFLICT or sub-locus
  if(s=="SEED_CONFLICT" || index(l,".")>0){
    wf=""; wv=0; sf2=""; sv2=0
    n=asorti(bits,ks)
    for(ki=1;ki<=n;ki++){
      c=ks[ki]; bv=bits[c]
      if(bv>wv){sf2=wf;sv2=wv;wf=c;wv=bv} else if(bv>sv2){sf2=c;sv2=bv}
    }
    base=(index(l,".")>0)?"SEGMENT_RESOLVED":"SEED_CONFLICT_RESOLVED"
    flag=(sv2>0 && (wv-sv2)<DA)?base"_AMBIGUOUS":base
    print l,wf,wv,sf2,sv2,(wv-sv2),flag; next
  }

  # Branch 2: NO_SEED
  if(s==""){
    wf=""; wv=0; sf2=""; sv2=0
    n=asorti(bits,ks)
    for(ki=1;ki<=n;ki++){
      c=ks[ki]; bv=bits[c]
      if(bv>wv){sf2=wf;sv2=wv;wf=c;wv=bv} else if(bv>sv2){sf2=c;sv2=bv}
    }
    flag=(sv2>0 && (wv-sv2)<DA)?"NO_SEED_AMBIGUOUS":"NO_SEED_RESOLVED"
    print l,wf,wv,sf2,sv2,(wv-sv2),flag; next
  }

  # Branch 3: UNIQUE_SEED
  sb=bits[s]; best=""; bestv=0
  for(c in bits) if(c!=s && bits[c]>bestv){best=c;bestv=bits[c]}
  delta=bestv-sb
  flag="UNIQUE_SEED_ONLY"
  if(bestv>0) flag="UNIQUE_SEED_WITH_COMPETITOR"
  if(delta>=D) flag="UNIQUE_SEED_CONTRADICTION"
  print l,s,sb,best,bestv,delta,flag

}' "$OUT/MASTER.locus_seeds.tsv" "$OUT/MASTER.locus_bits.tsv" \
> "$OUT/MASTER.final_assignment.raw.tsv"
check_file "$OUT/MASTER.final_assignment.raw.tsv"

awk 'BEGIN{FS=OFS="\t"}
FNR==NR{ prov[$1]=$2; next }
{ print $0, (($1 in prov)?prov[$1]:".") }
' "$OUT/MASTER.locus_provenance.segmented.tsv" \
  "$OUT/MASTER.final_assignment.raw.tsv" \
> "$OUT/MASTER.final_assignment.tsv"
check_file "$OUT/MASTER.final_assignment.tsv"

##############################################################################
# REPORT
##############################################################################
log "DONE"
echo ""
echo "=== Summary ==="
echo "Loci (original):  $(wc -l < "$OUT/MASTER.loci.bed")"
echo "Loci (segmented): $(wc -l < "$LOCI_BED")"
echo ""
echo "Family (col2):"; awk -F'\t' '{print $2}' "$OUT/MASTER.final_assignment.tsv" | sort | uniq -c | sort -rn
echo ""; echo "Flags (col7):"
awk -F'\t' '{print $7}' "$OUT/MASTER.final_assignment.tsv" | sort | uniq -c | sort -rn
echo ""
echo "Provenance (col8):"
awk -F'\t' '{print $8}' "$OUT/MASTER.final_assignment.tsv" | sort | uniq -c | sort -rn
echo ""
echo "Schema: locus_id|family|assigned_bits|best_other|best_other_bits|delta|flag|provenance"
echo "DELTA=$DELTA  DELTA_AMBIG=$DELTA_AMBIG"