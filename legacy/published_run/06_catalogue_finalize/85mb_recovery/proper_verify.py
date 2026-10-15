# -*- coding: utf-8 -*-
"""Apply the ESTABLISHED methodology (from chunk18_clsat_search.py, May 2026) to the region
this session excluded as 'telomeric' -- mask only the literal dense simple-repeat runs
(TTAGGG-family AND TAA-family, all rotations, both strands) and actively search the
REMAINING sequence for genuine CLsat (period 80-2000, matching this session's scan
parameters) rather than blanket-excluding the whole span."""
import re, os, subprocess

CHUNK44_FA = "/staging/tmp/clsat_v6/scan85/chunks/chunk_44.fa"
OUT_DIR = "/staging/tmp/clsat_v6/scan85/chunk44_properverify"
TRF = "/staging/miniconda3/bin/trf"
SCF = "JAWDKN010000003.1"
CHUNK_START = 83850001  # 1-based genomic start of chunk_44 (from chunk_coords.tsv)

MASK_MINLEN = 100
TRF_MAXP = 2000   # match this session's actual scan parameters, not the narrower 300
TRF_TIMEOUT = 900  # 15 min ceiling -- if masking truly worked this should finish fast

UNITS = ['TTAGGG', 'TAGGGT', 'AGGGTT', 'GGGTTA', 'GGTTAG', 'GTTAGG',
         'CCCTAA', 'CCTAAC', 'CTAACC', 'TAACCC', 'AACCCT', 'ACCCTA',
         'TAA', 'ATA', 'AAT', 'TTA', 'TAT', 'ATT']

os.makedirs(OUT_DIR, exist_ok=True)

with open(CHUNK44_FA) as fh:
    lines = fh.readlines()
seq = ''.join(l.strip() for l in lines[1:]).upper()
seqlen = len(seq)
print(f"chunk_44: {seqlen:,} bp (genomic {CHUNK_START:,}-{CHUNK_START+seqlen-1:,})")


def find_runs(seq, minlen):
    raw = []
    for unit in UNITS:
        pat = re.compile(f'(?:{unit}){{{max(2, minlen // len(unit))},}}', re.IGNORECASE)
        for m in pat.finditer(seq):
            if len(m.group()) >= minlen:
                raw.append((m.start(), m.end()))
    if not raw:
        return []
    raw.sort()
    merged = [raw[0]]
    for s, e in raw[1:]:
        if s <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))
    return merged


mask_runs = find_runs(seq, MASK_MINLEN)
total_masked = sum(e - s for s, e in mask_runs)
print(f"{len(mask_runs)} dense simple-repeat runs found, {total_masked:,} bp masked "
      f"({100*total_masked/seqlen:.1f}% of chunk)")
for s, e in mask_runs:
    print(f"  masked: genomic {CHUNK_START+s}-{CHUNK_START+e} ({e-s} bp)")

masked = bytearray(seq.encode())
for s, e in mask_runs:
    for i in range(s, e):
        masked[i] = ord('N')
masked = masked.decode()

mask_fa = os.path.join(OUT_DIR, "chunk44_masked.fa")
with open(mask_fa, 'w') as out:
    out.write(f">chunk44_masked\n")
    for i in range(0, len(masked), 80):
        out.write(masked[i:i+80] + '\n')

print(f"\nRunning TRF (MaxPeriod={TRF_MAXP}, timeout={TRF_TIMEOUT}s) on masked sequence...")
dat = f"{mask_fa}.2.5.7.80.10.50.{TRF_MAXP}.dat"
log = mask_fa + ".trf.log"
cmd = f"cd {OUT_DIR} && timeout {TRF_TIMEOUT} {TRF} {mask_fa} 2 5 7 80 10 50 {TRF_MAXP} -h -ngs > {dat} 2> {log}"
r = subprocess.run(cmd, shell=True)

if os.path.getsize(dat) == 0:
    print("TRF produced no output / may have timed out even after masking -- INVESTIGATE FURTHER")
else:
    hits = []
    with open(dat) as fh:
        for line in fh:
            if not line or not line[0].isdigit():
                continue
            p = line.split()
            if len(p) < 14:
                continue
            s, e, period, copies = int(p[0]), int(p[1]), int(p[2]), float(p[3])
            if 80 <= period <= 300 and copies >= 10.0:
                hits.append((CHUNK_START - 1 + s, CHUNK_START - 1 + e, period, copies, p[-2]))
    print(f"\nCLsat-range hits (period 80-300, copies>=10) in unmasked sequence: {len(hits)}")
    for gs, ge, period, copies, motif in sorted(hits):
        print(f"  {SCF}:{gs}-{ge}  period={period}  copies={copies:.1f}  motif={motif[:25]}")
    if not hits:
        print("CONFIRMED: no CLsat-range signal in chunk_44 outside the masked dense-repeat runs.")
