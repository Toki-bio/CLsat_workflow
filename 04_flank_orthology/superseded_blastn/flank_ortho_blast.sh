#!/bin/bash
set -euo pipefail
source /staging/conda/etc/profile.d/conda.sh 2>/dev/null || true
conda activate bioinfo 2>/dev/null || export PATH=/staging/conda/envs/bioinfo/bin:$PATH

DIR=/staging/Sergei/Lacertidae/Darevskia/Satellites/manuscript/CLsat_full_clusters_6genome
WORK=$DIR/flank_analysis
cd $WORK

# Split all_regions.fa into LF and RF sets, renaming headers to locus-centric IDs.
# We need: for each LF/RF region, which locus it belongs to. Use manifest.
# manifest: locus kind sp scf rs re trunc explen ; header = sp|scf:rs-re
awk -F'\t' '$2=="LF" && $5!=0{print $1"\t"$3"\t"$3"|"$4":"$5"-"$6}' manifest.tsv > lf_map.tsv
awk -F'\t' '$2=="RF" && $5!=0{print $1"\t"$3"\t"$3"|"$4":"$5"-"$6}' manifest.tsv > rf_map.tsv
# columns: locus  species  header

# Build python helper to emit per-side FASTA with headers "locus|species"
python3 - <<'PYEOF'
def load_fa(path):
    d={}; h=None; buf=[]
    for line in open(path):
        line=line.rstrip("\n")
        if line.startswith(">"):
            if h is not None: d[h]="".join(buf)
            h=line[1:]; buf=[]
        else: buf.append(line)
    if h is not None: d[h]="".join(buf)
    return d
seqs=load_fa("all_regions.fa")
for side in ("lf","rf"):
    out=open(f"{side}.fa","w")
    for line in open(f"{side}_map.tsv"):
        locus,sp,hdr=line.rstrip("\n").split("\t")
        s=seqs.get(hdr,"")
        if len(s)>=100:   # skip very short
            out.write(f">{locus}|{sp}\n{s}\n")
    out.close()
PYEOF

echo "LF records: $(grep -c '^>' lf.fa)"
echo "RF records: $(grep -c '^>' rf.fa)"

# All-vs-all blastn within each side
for side in lf rf; do
  makeblastdb -in ${side}.fa -dbtype nucl -out ${side}_db >/dev/null 2>&1
  blastn -query ${side}.fa -db ${side}_db \
    -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
    -evalue 1e-5 -max_target_seqs 50 -num_threads 8 \
    > ${side}.blast.tsv 2>/dev/null
  echo "${side} blast hits: $(wc -l < ${side}.blast.tsv)"
done

echo "DONE blast"
