"""Flank classification shared by clean_ortho_pairs.py and rebuild_ortho_segmented.py.

Classes, tested in this order:
  truncated   shorter than min_aln
  telomeric   >= 8 copies of the vertebrate telomere motif (TTAGGG or its reverse complement CCCTAA), or
              3-7 copies in a flank that is also low in complexity (a worn-down telomere stays a telomere)
  lowcomplex  simple sequence that is not telomeric: the share of the possible distinct 4-letter words that occur
              is below 0.5. "Possible" is min(256, len - 3): a 600-letter flank can hold at most 256 distinct words.
  unique      everything else

Fixed 2026-10-17. The earlier test divided the distinct-word count by the flank length (len - 3) and rejected below
0.35. In a 600-letter flank that ratio cannot exceed 256/597 = 0.43, so it rejected 65% of random 600-letter windows
of the D. unisexualis unm assembly (1 307 of 2 000). The scaled test rejects 2.5% (49 of 2 000).
"""

TELO_MOTIFS = ('TTAGGG', 'CCCTAA')
TELO_MIN = 8            # copies that make a flank telomeric on their own
TELO_WEAK = 3           # copies that make a simple flank telomeric
LOWC_MAX_SHARE = 0.5    # distinct 4-letter words / possible distinct words


def telo_copies(seq):
    s = seq.upper()
    return max(s.count(m) for m in TELO_MOTIFS)


def word_share(seq, k=4):
    s = seq.upper()
    n = len(s) - k + 1
    if n <= 0:
        return 0.0
    return len(set(s[i:i + k] for i in range(n))) / min(4 ** k, n)


def classify(seq, min_aln):
    if len(seq) < min_aln:
        return 'truncated'
    t = telo_copies(seq)
    simple = word_share(seq) < LOWC_MAX_SHARE
    if t >= TELO_MIN or (simple and t >= TELO_WEAK):
        return 'telomeric'
    if simple:
        return 'lowcomplex'
    return 'unique'
