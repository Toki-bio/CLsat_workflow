# -*- coding: utf-8 -*-
"""Classify the 85Mb-rescan TRF candidates (period in CLsat range, copies>=10) against the
four CLsat family consensus sequences by best rotation-aware alignment identity. Keep only
genuine matches; discard coincidental period matches to unrelated repeat families."""
import subprocess

DAT = '/staging/tmp/clsat_v6/scan85/scaf.trf.dedup.dat'
SD = '/staging/Sergei/Lacertidae/Darevskia/Satellites/motifs_output/clean/cross_clusters/selected/subclusters/seeds'

cons = {}
for fam in ('CLsat1', 'CLsat2', 'CLsat3', 'CLsat4'):
    seq = ''
    for line in open(f'{SD}/{fam}.cons.fa'):
        if not line.startswith('>'):
            seq += line.strip()
    cons[fam] = seq.upper()

def best_rotation_identity(motif, consensus):
    """Best identity of motif against any rotation of consensus (and its revcomp), same-length
    sliding compare -- handles TRF's arbitrary phase and strand for the reported motif."""
    def revcomp(s):
        comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
        return ''.join(comp.get(c, 'N') for c in reversed(s))

    best = 0.0
    for base in (consensus, revcomp(consensus)):
        doubled = base + base
        L = len(base)
        m = motif[:L] if len(motif) >= L else motif + motif[:L - len(motif)]
        for offset in range(L):
            window = doubled[offset:offset + L]
            if len(window) != len(m):
                continue
            matches = sum(1 for a, b in zip(window, m) if a == b)
            ident = matches / L
            if ident > best:
                best = ident
    return best

# candidates: period in CLsat range (114-118 or 143-148), copies>=10
candidates = []
for line in open(DAT):
    p = line.split()
    if len(p) < 14:
        continue
    s, e, period, copies = int(p[0]), int(p[1]), int(p[2]), float(p[3])
    motif = p[13].upper()
    if not (114 <= period <= 118 or 143 <= period <= 148):
        continue
    if copies < 10:
        continue
    candidates.append((s, e, period, copies, motif))

print(f'{len(candidates)} candidates with period in CLsat range and copies>=10\n')
kept = []
for s, e, period, copies, motif in candidates:
    scores = {fam: best_rotation_identity(motif, c) for fam, c in cons.items()}
    best_fam = max(scores, key=scores.get)
    best_score = scores[best_fam]
    verdict = 'KEEP' if best_score >= 0.70 else 'discard'
    print(f'{s}-{e}  period={period}  copies={copies:.1f}  best={best_fam} ({best_score:.2f})  {verdict}')
    if verdict == 'KEEP':
        kept.append((s, e, period, copies, best_fam, best_score))

print(f'\n{len(kept)} of {len(candidates)} candidates kept as genuine CLsat matches (>=70% identity)')
with open('/staging/tmp/clsat_v6/scan85/candidates_85mb.tsv', 'w') as out:
    out.write('start\tend\tperiod\tcopies\tfamily\tidentity\n')
    for s, e, period, copies, fam, score in kept:
        out.write(f'{s}\t{e}\t{period}\t{copies}\t{fam}\t{score:.3f}\n')
print('wrote candidates_85mb.tsv')
