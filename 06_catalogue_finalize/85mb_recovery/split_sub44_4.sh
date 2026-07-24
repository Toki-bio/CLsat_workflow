#!/usr/bin/env bash
set -euo pipefail
cd /staging/tmp/clsat_v6/scan85/chunks/sub44
samtools faidx sub44_4.fa

mkdir -p sub44_4_split
# sub44_4.fa length: 482922. Split into 4 pieces ~130kb with 10kb overlap.
declare -a REGIONS=("1-130000" "120001-250000" "240001-370000" "360001-482922")
declare -a OFFSETS=(0 120000 240000 360000)

for i in 0 1 2 3; do
  n=$((i+1))
  samtools faidx sub44_4.fa "sub44_4:${REGIONS[$i]}" | sed "1s/.*/>sub44_4_${n}/" > sub44_4_split/sub44_4_${n}.fa
done

echo "=== running TRF on 4 sub-sub-chunks (parallel) ==="
cd sub44_4_split
ls sub44_4_*.fa | parallel --jobs 4 'trf {} 2 5 7 80 10 50 2000 -h -ngs > {/.}.log 2>&1'
cd ..

echo "=== translating sub-sub coords back to sub44_4-local coords ==="
python3 - <<'PYEOF'
offsets = [0, 120000, 240000, 360000]
out_lines = []
for i in range(4):
    n = i + 1
    off = offsets[i]
    path = f"sub44_4_split/sub44_4_{n}.fa.2.5.7.80.10.50.2000.dat"
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

with open("sub44_4.fa.2.5.7.80.10.50.2000.dat", "w") as out:
    out.write("@sub44_4\n")
    out.writelines(dedup)

print(f"wrote sub44_4.fa.2.5.7.80.10.50.2000.dat: {len(dedup)} deduped rows from {len(out_lines)} raw rows")
PYEOF

touch sub44_4.split.done
echo DONE_SPLIT_SUB44_4
