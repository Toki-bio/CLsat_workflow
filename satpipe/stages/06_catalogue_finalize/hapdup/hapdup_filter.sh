#!/bin/bash
# Haplotype-duplicate rule (stage 16_regions, optional; on only when the config gives hapdup.pairs).
#
# Why: in a trio-phased assembly (hifiasm trio binning) reads that carry no parent-specific k-mers are
# placed in BOTH haplotypes. A pure satellite array has few such k-mers, so it can be assembled twice:
# once on a contig anchored by flanking sequence in its true haplotype, and once as a floating
# array-only contig in the other haplotype (found in D. unisexualis 224, 2026-11-08).
#
# Rule, applied in both directions of each pair a:b (a and b = the two haplotype copies of one animal):
#   (1) a contig of genome g is "array-only" when the catalogue loci on it cover all of it except at
#       most EDGE_TOL bp at each end;
#   (2) the array-only contig aligns to the partner genome p (minimap2 -x asm5 -c) over >= MIN_QCOV of
#       its length with gap-compressed identity (1 - de) >= MIN_IDENT (gap-compressed, so a shifted
#       repeat unit counts as one event, not as hundreds of letters);
#   (3) >= MIN_TOVL of the aligned partner stretch lies inside partner loci on contigs that are NOT
#       array-only (the twin is anchored).
#   -> every locus on that contig is REMOVED from locus_master.tsv.
# An array-only contig whose twin sits on an array-only partner contig is kept and flagged
# FLOATING_TWIN (neither copy is anchored, so the rule cannot tell which is the duplicate).
#
# Usage: hapdup_filter.sh <locus_master.tsv> <genomes.tsv> <out_dir> <pairs: a:b[,c:d...]>
# Env:   EDGE_TOL (100) MIN_QCOV (0.95) MIN_IDENT (0.99) MIN_TOVL (0.5) THREADS (8)
# Out:   <out_dir>/hapdup.tsv  one row per locus on an array-only contig, with evidence and verdict;
#        <locus_master.tsv> rewritten without REMOVED loci; the input kept as locus_master.prehapdup.tsv.
set -euo pipefail
export LC_ALL=C
[[ $# -eq 4 ]] || { echo "Usage: $0 <locus_master.tsv> <genomes.tsv> <out_dir> <pairs a:b[,c:d]>" >&2; exit 1; }
LM=$(readlink -f "$1"); GT=$(readlink -f "$2"); mkdir -p "$3"; OUT=$(readlink -f "$3"); PAIRS=$4
EDGE_TOL=${EDGE_TOL:-100}; MIN_QCOV=${MIN_QCOV:-0.95}; MIN_IDENT=${MIN_IDENT:-0.99}; MIN_TOVL=${MIN_TOVL:-0.5}
T=${THREADS:-8}
MM=${SATPIPE_TOOL_MINIMAP2:-minimap2}; ST=${SATPIPE_TOOL_SAMTOOLS:-samtools}
command -v "$MM" >/dev/null || { echo "FATAL: minimap2 not found (needed by hapdup)" >&2; exit 1; }
command -v "$ST" >/dev/null || { echo "FATAL: samtools not found" >&2; exit 1; }
W=$OUT/hapdup_work; mkdir -p "$W"; cd "$W" || exit 1

gpath() { awk -F'\t' -v g="$1" '$1==g{print $2; exit}' "$GT"; }

# contig -> "array-only" flag for genome $1 (needs its .fai): loci cover >= len - 2*EDGE_TOL
contig_flags() {
  awk -F'\t' -v g="$1" -v tol="$EDGE_TOL" 'FNR==NR{len[$1]=$2; next}
    $5==g { if(!($2 in len)){ print "FATAL: contig " $2 " of " g " not in fai" > "/dev/stderr"; bad=1; exit }
            cov[$2]+=$4-$3 }
    END{ if(bad) exit 1; for(c in cov) print c "\t" len[c] "\t" ((cov[c] >= len[c]-2*tol) ? 1 : 0) }' "$2.fai" "$LM" | sort
}

HDR='locus\tgenome\tcontig\tcontig_len\tpartner\tpartner_hit\tqcov\tident_gc\tpartner_loci\tverdict'
printf "$HDR\n" > "$OUT/hapdup.tsv"
IFS=, read -ra PP <<< "$PAIRS"
for pr in "${PP[@]}"; do
  [[ $pr == *:* ]] || { echo "FATAL: pair '$pr' is not a:b" >&2; exit 1; }
  for dir in "${pr%%:*} ${pr##*:}" "${pr##*:} ${pr%%:*}"; do
    read -r g p <<< "$dir"
    G=$(gpath "$g"); P=$(gpath "$p")
    [[ -n $G && -n $P ]] || { echo "FATAL: genome $g or $p not in $GT" >&2; exit 1; }
    for f in "$G" "$P"; do [[ -s $f.fai ]] || "$ST" faidx "$f"; done
    contig_flags "$g" "$G" > "$g.contigs.tsv"
    contig_flags "$p" "$P" > "$p.contigs.tsv"
    awk -F'\t' '$3==1{print $1}' "$g.contigs.tsv" > "$g.arrayonly.txt"
    if [[ ! -s $g.arrayonly.txt ]]; then echo "[hapdup] $g: no array-only contig"; continue; fi
    xargs "$ST" faidx "$G" < "$g.arrayonly.txt" > "$g.arrayonly.fa"
    "$MM" -c -x asm5 -N 50 -p 0.3 -t "$T" "$P" "$g.arrayonly.fa" > "${g}_vs_$p.paf" 2> "${g}_vs_$p.mm.log"
    # verdict per array-only contig, then one row per locus on it
    awk -F'\t' -v g="$g" -v p="$p" -v qc="$MIN_QCOV" -v mi="$MIN_IDENT" -v mt="$MIN_TOVL" '
      FILENAME==ARGV[1] { fl[$1]=$3; next }                                   # partner contig flags
      FILENAME==ARGV[2] { if($5==p){ n++; pc[n]=$2; ps[n]=$3; pe[n]=$4; pid[n]=$1 }
                          if($5==g){ loci[$2]=loci[$2] (loci[$2]?",":"") $1 }; next }
      FILENAME==ARGV[3] { clen[$1]=$2; if($3==1) ao[$1]=1; next }             # own contig flags
      {  q=$1; if(!(q in ao)) next
         de=-1; for(i=13;i<=NF;i++) if($i ~ /^de:f:/) de=substr($i,6)+0
         if(de<0) next
         cov=($4-$3)/$2; id=1-de; tl=$9-$8
         if(cov<qc || id<mi) {                     # failing hit: kept only as evidence (best one)
            if(rank[q]==0 && cov*id>bq[q]*bi[q]) { bq[q]=cov; bi[q]=id; bh[q]=$6":"$8"-"$9 }
            next }
         oa=0; of=0; ids=""
         for(k=1;k<=n;k++) if(pc[k]==$6){ s=(ps[k]>$8?ps[k]:$8); e=(pe[k]<$9?pe[k]:$9)
             if(e>s){ if(fl[$6]==1) of+=e-s; else oa+=e-s; ids=ids (ids?",":"") pid[k] } }
         v = (oa/tl>=mt) ? "REMOVED" : ((of/tl>=mt) ? "FLOATING_TWIN" : "TWIN_OUTSIDE_LOCI")
         r = (v=="REMOVED")?3:((v=="FLOATING_TWIN")?2:1)
         if(r>rank[q] || (r==rank[q] && cov*id>bq[q]*bi[q])) {
            rank[q]=r; best[q]=v; bq[q]=cov; bi[q]=id; bh[q]=$6":"$8"-"$9; bl[q]=ids }
      }
      END{ for(q in ao){ v=(best[q]!="")?best[q]:"NO_TWIN"
             split(loci[q], L, ",")
             for(j in L) printf "%s\t%s\t%s\t%d\t%s\t%s\t%.4f\t%.4f\t%s\t%s\n", L[j], g, q, clen[q], p,
                    ((q in bh)?bh[q]:"-"), ((q in bq)?bq[q]:0), ((q in bi)?bi[q]:0), (bl[q]?bl[q]:"-"), v } }
    ' "$p.contigs.tsv" "$LM" "$g.contigs.tsv" "${g}_vs_$p.paf" | sort -k1,1V >> "$OUT/hapdup.tsv"
    echo "[hapdup] $g vs $p: $(wc -l < "$g.arrayonly.txt") array-only contigs aligned"
  done
done

cp "$LM" "$(dirname "$LM")/locus_master.prehapdup.tsv"
awk -F'\t' 'FNR==NR{ if($10=="REMOVED") rm[$1]=1; next } !($1 in rm)' "$OUT/hapdup.tsv" "$(dirname "$LM")/locus_master.prehapdup.tsv" > "$LM.tmp"
mv "$LM.tmp" "$LM"
echo "[hapdup] verdicts:"; awk -F'\t' 'NR>1{print $2"\t"$10}' "$OUT/hapdup.tsv" | sort | uniq -c
echo "[hapdup] locus_master: $(wc -l < "$(dirname "$LM")/locus_master.prehapdup.tsv") -> $(wc -l < "$LM") loci"
