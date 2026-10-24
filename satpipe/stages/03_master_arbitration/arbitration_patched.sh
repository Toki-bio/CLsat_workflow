#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C   # byte-order sort: locus numbering must not depend on the locale

# arbitration_patched.sh -- Step 14: locus segmentation and family assignment.
#
# Usage: arbitration_patched.sh <families.tsv> <MASTER.merged.bed> <out_dir>
#   families.tsv       one row per family, tab-separated, '#' lines ignored:
#                        family  MONO  <family>.bits.tsv (stage-13 rescan of MASTER)  <family> seed BED (stage 10)
#                      The row order fixes the column order of MASTER.locus_bits.tsv.
#   MASTER.merged.bed  stage-12 output (chrom start end tag [provenance]); locus L<n> = row n; the tag
#                      starts with the genome code (`sp|scf:s-e|...`), as do the seed and bits rows
#   out_dir            every output goes here (created)
# Environment (defaults = published Darevskia run):
#   DELTA          UNIQUE_SEED contradiction threshold                       [20]
#   DELTA_AMBIG    CONFLICT / sub-locus / NO_SEED ambiguity threshold          [15]
#   SEG_BIN        D3b bin width, bp                                           [1000]
#   SEG_MIN_RUN    bins a competing family must dominate to open a segment     [10]
#   SEG_MIN_LOCUS  segments shorter than this are dropped or merged (SEG_SHORT) [2000]
#   SEG_SHORT      drop (published) | merge                                    [drop]
#   SCAN_WIN_MONO  stage-13 scan window in monomers (D3 width = this x MONO)  [2]
#
# D3 converts each scan row (1-based window start `pos`) to a 0-based BED interval of width
# 2 x MONO, i.e. the scan window of stage 13 (292/290/292/232 bp for the Darevskia families).
#
# Patches vs original:
#  D3 : pos-1 for 0-based BED start; window = 2 x MONO
#  D3b: sorted() tiebreak in Python segmentation (deterministic); bins on pos-1 like D3
#  D6 : three explicit branches: UNIQUE_SEED / SEED_CONFLICT+sub / NO_SEED
#       sub-loci and SEED_CONFLICT get delta scrutiny via DELTA_AMBIG
#       NO_SEED uses bitscore winner without phantom 0-delta comparison
#       branch 3 writes delta = seed - strongest other (same sign as branches 1-2)
#
# Output columns (MASTER.final_assignment.tsv): locus_id family assigned_bits best_other
#                 best_other_bits delta flag provenance(col8)
#
# Genomes may share contig names (two hifiasm assemblies both have h1tg000102l). Every interval is
# therefore keyed internally by `genome|contig` (the *.key.bed files, SEEDS.tagged.bed, <fam>.bits.bed),
# so seeds, scan windows and segmentation bins of one genome never reach a locus of another. The public
# MASTER.loci.bed and MASTER.loci.segmented.bed keep contig start end locus and add the genome as col 5
# (.short_segments.bed: col 7).

die(){ echo "FATAL: $*" >&2; exit 1; }
check_file(){ [[ -s "$1" ]] || die "Missing: $1"; }
log(){ echo "[$(date +"%F %T")] $*"; }

[[ $# -eq 3 ]] || die "Usage: $0 <families.tsv> <MASTER.merged.bed> <out_dir>"
FAMTAB="$1"; MASTER_BED="$2"; OUT="$3"
check_file "$FAMTAB"; check_file "$MASTER_BED"
mkdir -p "$OUT"

DELTA="${DELTA:-20}"              # UNIQUE_SEED contradiction threshold
DELTA_AMBIG="${DELTA_AMBIG:-15}"  # CONFLICT/sub-locus/NO_SEED ambiguity threshold
SEG_BIN="${SEG_BIN:-1000}"; SEG_MIN_RUN="${SEG_MIN_RUN:-10}"; SEG_MIN_LOCUS="${SEG_MIN_LOCUS:-2000}"
export SEG_SHORT="${SEG_SHORT:-drop}"
SCAN_WIN_MONO="${SCAN_WIN_MONO:-2}"   # stage-13 scan window in monomers (scan.win_mono); D3 window = this x MONO
[[ "$SCAN_WIN_MONO" =~ ^[1-9][0-9]*$ ]] || die "SCAN_WIN_MONO must be a positive integer, got '$SCAN_WIN_MONO'"

# families table -> parallel arrays
declare -a FAMS=() MONOS=() BITS=() SEEDS=()
while IFS=$'\t' read -r fam mono bits seed _; do
  [[ -z "$fam" || "$fam" == \#* ]] && continue
  [[ "$mono" =~ ^[0-9]+$ ]] || die "families.tsv: MONO must be an integer for $fam, got '$mono'"
  FAMS+=("$fam"); MONOS+=("$mono"); BITS+=("$bits"); SEEDS+=("$seed")
done < "$FAMTAB"
NF_=${#FAMS[@]}
(( NF_ >= 1 )) || die "families.tsv has no rows"
FAMS_STR="${FAMS[*]}"

# provenance column present?
if awk -F'\t' 'NR==1{exit !(NF>=5)}' "$MASTER_BED"; then HAS_PROV=1; else HAS_PROV=0; fi

log "Pre-flight (families=${FAMS_STR} DELTA=${DELTA} DELTA_AMBIG=${DELTA_AMBIG} SEG_SHORT=${SEG_SHORT} prov=${HAS_PROV})"
for i in "${!FAMS[@]}"; do check_file "${BITS[$i]}"; check_file "${SEEDS[$i]}"; done
command -v bedtools >/dev/null || die "bedtools not in PATH"
command -v python3  >/dev/null || die "python3 not in PATH"
command -v gawk     >/dev/null || die "gawk not in PATH (asorti)"

##############################################################################
# D1
##############################################################################
log "D1: assigning locus IDs"

# genome code = the tag up to its first '|'
awk 'BEGIN{FS=OFS="\t"} { i=index($4,"|"); if(i<2){ print "FATAL: MASTER row " NR " has no genome code in col 4: " $4 > "/dev/stderr"; exit 1 }
  sp=substr($4,1,i-1); print sp"|"$1,$2,$3,"L"NR > K; print $1,$2,$3,"L"NR,sp }' K="$OUT/MASTER.loci.key.bed" "$MASTER_BED" > "$OUT/MASTER.loci.bed"

if [[ "$HAS_PROV" == "1" ]]; then
  awk 'BEGIN{FS=OFS="\t"}{print "L"NR, ($5?$5:".")}' "$MASTER_BED" \
    > "$OUT/MASTER.locus_provenance.tsv"
else
  awk 'BEGIN{OFS="\t"}{print "L"NR, "."}' "$MASTER_BED" \
    > "$OUT/MASTER.locus_provenance.tsv"
fi
check_file "$OUT/MASTER.loci.bed"; check_file "$OUT/MASTER.loci.key.bed"

##############################################################################
# D2
##############################################################################
log "D2: tagging seeds"

for i in "${!FAMS[@]}"; do
  awk -v F="${FAMS[$i]}" 'BEGIN{FS=OFS="\t"}{ i=index($4,"|"); if(i<2){ print "FATAL: seed row without genome code: " $4 > "/dev/stderr"; exit 1 }
    print substr($4,1,i-1)"|"$1,$2,$3,F }' "${SEEDS[$i]}"
done | sort -k1,1 -k2,2n > "$OUT/SEEDS.tagged.bed"
check_file "$OUT/SEEDS.tagged.bed"

##############################################################################
# D3 -- PATCHED: 0-based BED start; window = 2 x MONO
##############################################################################
log "D3: bits TSV -> BED (pos-1; width 2 x MONO)"

for i in "${!FAMS[@]}"; do
  w=$(( SCAN_WIN_MONO * MONOS[i] ))   # must equal the stage-13 window (--scan-win-bp = scan.win_mono x MONO)
  awk -v W="$w" 'BEGIN{FS=OFS="\t"} NR>1{ i=index($4,"|"); if(i<2){ print "FATAL: bits row without genome code: " $4 > "/dev/stderr"; exit 1 }
    print substr($4,1,i-1)"|"$1, $2-1, $2-1+W, $3 }' "${BITS[$i]}" > "$OUT/${FAMS[$i]}.bits.bed"
  check_file "$OUT/${FAMS[$i]}.bits.bed"
done

##############################################################################
# D3b -- PATCHED: deterministic sorted() tiebreak
##############################################################################
log "D3b: segmenting composite loci"

FAMBITS_ARGS=()
for i in "${!FAMS[@]}"; do FAMBITS_ARGS+=("${FAMS[$i]}" "${BITS[$i]}"); done
python3 - "$OUT/MASTER.loci.key.bed" "$SEG_BIN" "$SEG_MIN_RUN" "$SEG_MIN_LOCUS" \
  "$OUT/MASTER.loci.segmented.key.bed" "${FAMBITS_ARGS[@]}" <<'PYEOF'
import sys, os
from collections import defaultdict
loci_file=sys.argv[1]
BIN=int(sys.argv[2]); MIN_RUN=int(sys.argv[3])
MIN_LOCUS=int(sys.argv[4]); out_file=sys.argv[5]
rest=sys.argv[6:]
families=rest[0::2]; bits_files=rest[1::2]
loci=[]
with open(loci_file) as f:
    for line in f:
        p=line.rstrip("\n").split("\t"); loci.append((p[0],int(p[1]),int(p[2]),p[3]))
print(f"  Loaded {len(loci)} loci",file=sys.stderr)
bins=defaultdict(lambda:defaultdict(lambda:defaultdict(float)))
for fam,bf in zip(families,bits_files):
    with open(bf) as f:
        f.readline()
        for line in f:
            p=line.rstrip("\n").split("\t")
            scf,pos,bits=p[3].split("|",1)[0]+"|"+p[0],int(p[1]),float(p[2])   # genome|contig, as the loci
            b=(pos-1)//BIN   # pos is a 1-based window start; bins are 0-based like the BED (same pos-1 as D3)
            if bits>bins[scf][b][fam]: bins[scf][b][fam]=bits
    print(f"  Loaded {bf} ({fam})",file=sys.stderr)
total_seg=0; n_split=0; out=open(out_file,"w")
# Segments shorter than MIN_LOCUS: the published run DROPPED them (default, SEG_SHORT=drop),
# which removes sequence from the catalogue and leaves gaps in segment numbering
# (e.g. L12.1, L12.3). SEG_SHORT=merge folds each short segment into its longer neighbour
# instead. Either way every dropped/merged piece is written to <out>.short_segments.bed.
SEG_SHORT=os.environ.get("SEG_SHORT","drop")
short_out=open(out_file+".short_segments.bed","w"); n_short=0; bp_short=0
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
        segs=[]
        for si,(sb,eb,fam) in enumerate(segments):
            s_bp=max(start,(bin_start+sb)*BIN); e_bp=min(end,(bin_start+eb+1)*BIN)
            segs.append([s_bp,e_bp,fam,si+1])
        if SEG_SHORT=="merge":
            # fold short segments into the longer adjacent segment (keep original numbering)
            changed=True
            while changed and len(segs)>1:
                changed=False
                for k,(s_bp,e_bp,fam,num) in enumerate(segs):
                    if e_bp-s_bp<MIN_LOCUS:
                        short_out.write(f"{scf}\t{s_bp}\t{e_bp}\t{lid}.{num}\t{fam}\tmerged\n"); n_short+=1; bp_short+=e_bp-s_bp
                        nb=k-1 if k>0 and (k==len(segs)-1 or (segs[k-1][1]-segs[k-1][0])>=(segs[k+1][1]-segs[k+1][0])) else k+1
                        segs[nb][0]=min(segs[nb][0],s_bp); segs[nb][1]=max(segs[nb][1],e_bp)
                        del segs[k]; changed=True; break
        if len(segs)==1:
            out.write(f"{scf}\t{segs[0][0]}\t{segs[0][1]}\t{lid}\n"); total_seg+=1; continue
        for s_bp,e_bp,fam,num in segs:
            if e_bp-s_bp<MIN_LOCUS:
                short_out.write(f"{scf}\t{s_bp}\t{e_bp}\t{lid}.{num}\t{fam}\tdropped\n"); n_short+=1; bp_short+=e_bp-s_bp
                continue
            out.write(f"{scf}\t{s_bp}\t{e_bp}\t{lid}.{num}\n"); total_seg+=1
out.close(); short_out.close()
print(f"  Original: {len(loci)}, split: {n_split}, output: {total_seg}",file=sys.stderr)
print(f"  Short segments (<{MIN_LOCUS} bp) {SEG_SHORT}: {n_short} pieces, {bp_short} bp -> {out_file}.short_segments.bed",file=sys.stderr)
PYEOF

check_file "$OUT/MASTER.loci.segmented.key.bed"
LOCI_BED="$OUT/MASTER.loci.segmented.key.bed"
# public tables: contig start end locus genome (short segments: contig start end locus family action genome)
awk 'BEGIN{FS=OFS="\t"}{ i=index($1,"|"); sp=substr($1,1,i-1); $1=substr($1,i+1); print $0,sp }' \
  "$LOCI_BED" > "$OUT/MASTER.loci.segmented.bed"
awk 'BEGIN{FS=OFS="\t"}{ i=index($1,"|"); sp=substr($1,1,i-1); $1=substr($1,i+1); print $0,sp }' \
  "$LOCI_BED.short_segments.bed" > "$OUT/MASTER.loci.segmented.bed.short_segments.bed"
rm -f "$LOCI_BED.short_segments.bed"
check_file "$OUT/MASTER.loci.segmented.bed"

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
log "  SEED_CONFLICT: $(grep -c SEED_CONFLICT "$OUT/MASTER.locus_seeds.tsv" || true)"

##############################################################################
# D4
##############################################################################
log "D4: max bits per locus"

for i in "${!FAMS[@]}"; do
  fam="${FAMS[$i]}"
  bedtools intersect -a "$LOCI_BED" -b "$OUT/${fam}.bits.bed" -wa -wb \
  | awk 'BEGIN{FS=OFS="\t"}
    { l=$4; b=$8; if(b>max[l]) max[l]=b }
    END{ for(l in max) print l,max[l] }' \
  > "$OUT/MASTER.${fam}.maxbits.tsv"
  check_file "$OUT/MASTER.${fam}.maxbits.tsv"
done

##############################################################################
# D5
##############################################################################
log "D5: locus bits matrix (columns: locus ${FAMS_STR})"

cut -f4 "$LOCI_BED" | awk -v N="$NF_" 'BEGIN{OFS="\t"}{s=$1; for(i=1;i<=N;i++) s=s OFS 0; print s}' > "$OUT/MASTER.locus_bits.tsv"

for i in "${!FAMS[@]}"; do
  col=$((i+2))
  awk -v C="$col" 'BEGIN{FS=OFS="\t"}
  FNR==NR{b[$1]=$2; next}
  { if($1 in b) $C=b[$1]; print }
  ' "$OUT/MASTER.${FAMS[$i]}.maxbits.tsv" "$OUT/MASTER.locus_bits.tsv" \
  > "$OUT/tmp" && mv "$OUT/tmp" "$OUT/MASTER.locus_bits.tsv"
done

##############################################################################
# D6 -- PATCHED: three explicit branches
##############################################################################
log "D6: arbitration (three-branch)"

gawk -v D="$DELTA" -v DA="$DELTA_AMBIG" -v FAMS="$FAMS_STR" 'BEGIN{FS=OFS="\t"; nf=split(FAMS,fam," ")}

FNR==NR{ if(FNR==1) next; seed[$1]=$2; next }

{
  l=$1
  delete bits
  for(fi=1;fi<=nf;fi++) bits[fam[fi]]=$(fi+1)+0
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
  # delta = seed bits - strongest other family, the same sign convention as branches 1-2
  # (assigned - runner-up). Before 2026-10-15 this branch wrote bestv-sb, so the column had
  # mixed signs. Candidates are visited in sorted order so ties are deterministic.
  sb=bits[s]; best=""; bestv=0
  n=asorti(bits,ks)
  for(ki=1;ki<=n;ki++){ c=ks[ki]; if(c!=s && bits[c]>bestv){best=c;bestv=bits[c]} }
  delta=sb-bestv
  flag="UNIQUE_SEED_ONLY"
  if(bestv>0) flag="UNIQUE_SEED_WITH_COMPETITOR"
  if(delta<=-D) flag="UNIQUE_SEED_CONTRADICTION"   # the seed family is retained; the flag marks the disagreement
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
echo "DELTA=$DELTA  DELTA_AMBIG=$DELTA_AMBIG  SEG_SHORT=$SEG_SHORT  families=$FAMS_STR"
