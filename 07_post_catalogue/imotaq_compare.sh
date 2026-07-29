#!/usr/bin/env bash
# =============================================================================
# Is DarSat2 the previously described IMO-TaqI satellite?
#
# The manuscript originally called this a candidate identity "pending direct
# sequence confirmation". This script is that confirmation.
#
# METHOD
#   1. fetch every GenBank record matching IMO-TaqI (220 records, 10 species);
#   2. compare them against all 39 DarSat family consensuses with ssearch36,
#      subjects doubled head-to-tail so the comparison is insensitive to monomer
#      phase;
#   3. take the best-scoring DarSat family per clone and tally.
#
# RESULT AS RUN (2026-07-28)
#   All 220 clones take DarSat2 as their top hit, up to 90.8% identity over 192 bp
#   (189.4 bits). Next best of the 39 families is DarSat34 at 35.0 bits; DarSat1 is
#   noise (24.2 bits over 55 bp). Per species, best identity 82-91% across Lacerta,
#   Timon, Gallotia, Atlantolacerta and Trogonophis.
#
#   => DarSat2 IS IMO-TaqI. This also overturned a note in the project's own
#      RESULTS_NUMBERS_LOCKED.md claiming the identity belonged to DarSat1.
#
# THE GALLOTIA CHECK (step 4, optional)
#   IMO-TaqI clones exist from Gallotia galloti, yet the outgroup screen scored 0
#   DarSat2 loci there. Both routes were checked and agree that this is a true
#   negative for that assembly, not a detection failure:
#     - TRF motif set (416 motifs): best DarSat2 match 27.2 bits over 62 bp,
#       against a family threshold of T(DarSat2) = 106.5;
#       Podarcis muralis positive control gives 131.1 bits over 177 bp.
#     - whole assembly GCA_042478105.1 (1.19 Gb): 4 hits, longest 65 bp, best
#       83.3% over 54 bp; nothing >=100 bp, nothing >=85%.
#   Consistent with pericentromeric heterochromatin being under-represented in a
#   contig-only assembly, or the clones coming from a different population.
# =============================================================================
set -euo pipefail

SS=/staging/conda/envs/bioinfo/bin/ssearch36
DARSAT=DarSat_all_consensus.fa          # the 39 family consensuses (Additional file 2)
OUT=${1:-./imotaq}
E=https://eutils.ncbi.nlm.nih.gov/entrez/eutils
mkdir -p "$OUT"; cd "$OUT"

# --- 1. fetch IMO-TaqI from GenBank -----------------------------------------
curl -s "$E/esearch.fcgi?db=nuccore&term=IMO-TaqI&retmax=300&usehistory=y" > es.xml
KEY=$(grep -o '<QueryKey>[0-9]*' es.xml | tr -d '<QueryKey>')
WEB=$(grep -o '<WebEnv>[^<]*' es.xml | sed 's/<WebEnv>//')
curl -s "$E/efetch.fcgi?db=nuccore&query_key=$KEY&WebEnv=$WEB&rettype=fasta&retmode=text" \
     > imotaq.fasta
echo "IMO-TaqI records fetched: $(grep -c '^>' imotaq.fasta)"
echo "species represented:"
grep '^>' imotaq.fasta | sed 's/^>[^ ]* //; s/ clone.*//' | sort | uniq -c | sort -rn

# --- 2. double the DarSat consensuses (phase-insensitive subjects) -----------
awk '/^>/{if(s)print h"\n"s s; h=$0; s=""; next}{s=s$0}END{print h"\n"s s}' \
    "$DARSAT" > darsat_doubled.fa

# --- 3. all-vs-all and tally the best family per clone ----------------------
"$SS" -n -m 8C -E 10 -b 5000 imotaq.fasta darsat_doubled.fa 2>/dev/null > hits.tsv

echo
echo "best DarSat family per IMO-TaqI clone, tallied:"
awk '!/^#/{if($12>b[$1]){b[$1]=$12; f[$1]=$2; p[$1]=$3}}
     END{for(q in b) printf "%s\t%.1f\t%.1f\n", f[q], b[q], p[q]}' hits.tsv \
 | awk '{n[$1]++; if($2>mx[$1])mx[$1]=$2; if($3>mp[$1])mp[$1]=$3}
        END{for(k in n) printf "  %-14s %3d clones   max bits %6.1f   max identity %.1f%%\n",
                                k, n[k], mx[k], mp[k]}' | sort -k2,2nr

echo
echo "best alignment per DarSat family (any clone):"
awk '!/^#/{if($12>b[$2]){b[$2]=$12; p[$2]=$3; a[$2]=$4}}
     END{for(s in b) printf "  %-14s bits %6.1f  identity %5.1f%%  aln %3d\n", s, b[s], p[s], a[s]}' \
    hits.tsv | sort -k3,3nr | head -6

# --- 4. optional: is the family really absent from a given assembly? --------
# usage: bash imotaq_compare.sh <outdir> <genome.fna> <DarSat2.fa>
if [ "${2:-}" ] && [ "${3:-}" ]; then
  echo
  echo "whole-genome check of $(basename "$3") against $(basename "$2"):"
  makeblastdb -in "$2" -dbtype nucl -out gdb >/dev/null 2>&1
  blastn -task blastn -word_size 11 -query "$3" -db gdb -evalue 1e-3 \
         -outfmt "6 qseqid sseqid pident length evalue bitscore" 2>/dev/null \
    | sort -k6,6gr | head
fi
