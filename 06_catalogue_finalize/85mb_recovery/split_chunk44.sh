#!/usr/bin/env bash
set -euo pipefail
cd /staging/tmp/clsat_v6/scan85/chunks
samtools faidx chunk_44.fa

mkdir -p sub44
# 4 sub-regions of chunk_44 (1-based), 50kb overlap
declare -a REGIONS=("1-550000" "500001-1050000" "1000001-1550000" "1500001-1982912")
declare -a OFFSETS=(0 500000 1000000 1500000)

for i in 0 1 2 3; do
  n=$((i+1))
  samtools faidx chunk_44.fa "chunk_44:${REGIONS[$i]}" | sed "1s/.*/>sub44_${n}/" > sub44/sub44_${n}.fa
done

echo "=== running TRF on 4 sub-chunks (parallel) ==="
cd sub44
ls sub44_*.fa | parallel --jobs 4 'trf {} 2 5 7 80 10 50 2000 -h -ngs > {/.}.log 2>&1'
cd ..

echo "=== translating sub-chunk-local coords back to chunk_44-local coords ==="
python3 - <<'PYEOF'
offsets = [0, 500000, 1000000, 1500000]
out_lines = []
for i in range(4):
    n = i + 1
    off = offsets[i]
    path = f"sub44/sub44_{n}.fa.2.5.7.80.10.50.2000.dat"
    try:
        f = open(path)
    except FileNotFoundError:
        print(f"WARNING: {path} not found, skipping")
        continue
    for line in f:
        if line.startswith('@'):
            continue
        p = line.split()
        if len(p) < 2:
            continue
        s, e = int(p[0]) + off, int(p[1]) + off
        newline = f"{s} {e} " + " ".join(p[2:]) + "\n"
        out_lines.append((s, e, newline))
    f.close()

out_lines.sort(key=lambda x: (x[0], x[1]))
dedup = []
prevkey = None
prevend = -10**9
for s, e, line in out_lines:
    key = (e - s)
    if key != prevkey or s - prevend > 100:
        dedup.append(line)
        prevkey = key
    prevend = e

with open("chunk_44.dat", "w") as out:
    out.write("@chunk_44\n")
    out.writelines(dedup)

print(f"wrote chunk_44.dat: {len(dedup)} deduped rows from {len(out_lines)} raw rows")
PYEOF

touch chunk_44.rerun.done
echo DONE_SPLIT44
