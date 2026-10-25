"""Anchor phasing: put every written consensus in one reading frame.

A circular monomer has no natural start. Step 5 rotates motifs to the lexicographically smallest
rotation (arbitrary phase, only for clustering); the manuscript figures and the deposited
consensuses use a common frame set by an anchor family (Methods 3.4: doubled consensus aligned to
the CLsat1 anchor, local alignment, start checked against the motif AAAGAAG). This module makes that
step automatic: `rotate_to_anchor` finds the rotation of a sequence whose start corresponds to the
anchor's first base, `check_frame` verifies it, and stage 21_anchor / `satpipe anchor` apply both to
every consensus the run writes, failing when a sequence cannot be placed or does not start at the
anchor motif.
"""
from __future__ import annotations

import os

MATCH, MISMATCH, GAP = 2, -3, -5          # local alignment scoring (nucleotides, linear gaps)


class AnchorError(Exception):
    pass


def read_fasta(path: str) -> dict:
    """Records by name; a name that appears twice is an error (it used to overwrite the first silently)."""
    out, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n\r')
            if line.startswith('>'):
                if name is not None:
                    out[name] = ''.join(buf)
                name, buf = line[1:].split()[0], []
                if name in out:
                    raise ValueError('%s: record name %s appears twice' % (path, name))
            elif name is not None:
                buf.append(line.strip())
    if name is not None:
        out[name] = ''.join(buf)
    return out


def write_fasta(path: str, recs: dict, width: int = 0) -> None:
    with open(path, 'w') as fh:
        for n, s in recs.items():
            fh.write('>%s\n' % n)
            if width:
                for i in range(0, len(s), width):
                    fh.write(s[i:i + width] + '\n')
            else:
                fh.write(s + '\n')


def clean(seq: str) -> str:
    return seq.replace('-', '').replace(' ', '').upper().replace('U', 'T')


def local_align(q: str, t: str):
    """Smith-Waterman, linear gaps. Returns (score, q_start, t_start, aligned_cols, matches); starts are
    0-based positions where the best local alignment begins."""
    n, m = len(q), len(t)
    prev = [0] * (m + 1)
    # trace[i][j]: 0 stop, 1 diag, 2 up (gap in t), 3 left (gap in q)
    trace = [bytearray(m + 1) for _ in range(n + 1)]
    best, bi, bj = 0, 0, 0
    for i in range(1, n + 1):
        cur = [0] * (m + 1)
        qi = q[i - 1]
        row = trace[i]
        for j in range(1, m + 1):
            d = prev[j - 1] + (MATCH if qi == t[j - 1] else MISMATCH)
            u = prev[j] + GAP
            l = cur[j - 1] + GAP
            s = max(0, d, u, l)
            cur[j] = s
            if s == 0:
                row[j] = 0
            elif s == d:
                row[j] = 1
            elif s == u:
                row[j] = 2
            else:
                row[j] = 3
            if s > best:
                best, bi, bj = s, i, j
        prev = cur
    i, j, cols, mat = bi, bj, 0, 0
    while i > 0 and j > 0 and trace[i][j]:
        k = trace[i][j]
        cols += 1
        if k == 1:
            mat += q[i - 1] == t[j - 1]
            i, j = i - 1, j - 1
        elif k == 2:
            i -= 1
        else:
            j -= 1
    return best, i, j, cols, mat


def hamming_prefix(seq: str, motif: str) -> int:
    return sum(1 for a, b in zip(seq, motif) if a != b) + max(0, len(motif) - len(seq))


def anchor_frame(anchor_seq: str, motif: str = '', max_mm: int = 1) -> str:
    """The anchor family's consensus, rotated so that it starts at `motif` (best-matching position in the
    doubled sequence; at most `max_mm` mismatches). Without a motif the sequence is used as given."""
    a = clean(anchor_seq)
    if not motif:
        return a
    motif = motif.upper()
    d = a + a
    cands = sorted((hamming_prefix(d[i:i + len(motif)], motif), i) for i in range(len(a)))
    mm, pos = cands[0]
    if mm > max_mm:
        raise AnchorError('anchor motif %s not found in the anchor consensus (best %d mismatches)' % (motif, mm))
    if len(cands) > 1 and cands[1][0] == mm:
        raise AnchorError('anchor motif %s occurs more than once in the anchor consensus at equal quality '
                          '(positions %d and %d): the frame is ambiguous' % (motif, pos, cands[1][1]))
    return d[pos:pos + len(a)]


def rotate_to_anchor(seq: str, anchor: str, min_identity: float = 0.6, min_aligned: int = 40,
                     motif: str = '', max_mm: int = 1, max_snap: int = 40):
    """Rotate `seq` so that it starts at the position corresponding to anchor[0].
    The alignment gives the start; when the anchor's first bases fall in a deleted or diverged stretch
    (e.g. CLsat4, 116 bp against a 146 bp anchor) that start is extrapolated and may be off, so with a
    `motif` the start is snapped to the nearest circular position that carries the motif (at most
    `max_snap` bases away, at most `max_mm` mismatches). Returns (rotated, shift, identity, aligned_cols,
    snapped_by) with snapped_by = 0 when the alignment start already carried the motif."""
    s = clean(seq)
    d = s + s
    score, qs, ts, cols, mat = local_align(anchor, d)
    ident = mat / cols if cols else 0.0
    if cols < min(min_aligned, len(anchor)) or ident < min_identity:
        raise AnchorError('cannot place the sequence on the anchor (aligned %d columns, %.0f%% identity)'
                          % (cols, 100 * ident))
    n = len(s)
    shift, snapped = (ts - qs) % n, 0
    if motif and hamming_prefix(d[shift:shift + len(motif)], motif.upper()) > max_mm:
        best = None      # exact motif matches first, then the nearest to the alignment start
        for p in range(n):
            mm = hamming_prefix(d[p:p + len(motif)], motif.upper())
            if mm <= max_mm:
                dist = min((p - shift) % n, (shift - p) % n)
                if dist <= max_snap and (best is None or (mm, dist) < best[:2]):
                    best = (mm, dist, p)
        if best is None:
            raise AnchorError('the motif %s is not within %d bases of the alignment start%s'
                              % (motif, max_snap, ''))   # no 'nearest' is computed (the message used to claim one, #118)
        shift, snapped = best[2], best[1]
    return d[shift:shift + n], shift, ident, cols, snapped


def check_frame(seq: str, motif: str, max_mm: int = 1) -> bool:
    if not motif:
        return True
    return hamming_prefix(clean(seq), motif.upper()) <= max_mm


def anchor_records(recs: dict, anchor: str, motif: str = '', max_mm: int = 1, min_identity: float = 0.6):
    """Rotate every record into the anchor frame. Returns (rotated dict, report rows). Raises AnchorError
    listing every record that could not be placed or does not start at the motif."""
    out, rows, bad = {}, [], []
    for name, seq in recs.items():
        try:
            rot, shift, ident, cols, snap = rotate_to_anchor(seq, anchor, min_identity=min_identity, motif=motif, max_mm=max_mm)
        except AnchorError as e:
            bad.append('%s: %s' % (name, e))
            rows.append((name, len(clean(seq)), '', '', '', 'unplaced'))
            continue
        ok = check_frame(rot, motif, max_mm)
        if not ok:
            bad.append('%s: after rotation the sequence starts %s, not the anchor motif %s' % (name, rot[:len(motif)], motif))
        out[name] = rot
        rows.append((name, len(rot), shift, '%.1f' % (100 * ident), rot[:max(7, len(motif))],
                     ('ok' if ok else 'motif-mismatch') + ('; snapped %d to motif' % snap if snap else '')))
    if bad:
        raise AnchorError('anchor phasing failed:\n  ' + '\n  '.join(bad))
    return out, rows


def verify_records(recs: dict, anchor: str, motif: str = '', max_mm: int = 1, min_identity: float = 0.6):
    """Check (without changing) that each record is already in the anchor frame (gaps are ignored).
    Returns a list of problem strings; empty means every record starts at the anchor."""
    problems = []
    for name, seq in recs.items():
        s = clean(seq)
        try:
            _rot, shift, _i, _c, _snap = rotate_to_anchor(s, anchor, min_identity=min_identity, motif=motif, max_mm=max_mm)
        except AnchorError as e:
            problems.append('%s: %s' % (name, e))
            continue
        if shift != 0 or not check_frame(s, motif, max_mm):
            problems.append('%s: out of the anchor frame (rotation needed: %d bases; starts %s)' % (name, shift, s[:max(7, len(motif))]))
    return problems


def fasta_files(path: str) -> list:
    """A FASTA file, or the FASTA files of a directory."""
    if os.path.isdir(path):
        # *.anchored.fa are this module's own outputs: re-reading them re-anchored stale copies (GLM audit #98)
        return sorted(os.path.join(path, f) for f in os.listdir(path)
                      if f.lower().endswith(('.fa', '.fas', '.fasta', '.fna')) and '.anchored.' not in f.lower())
    return [path]
