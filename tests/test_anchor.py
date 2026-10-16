"""Anchor phasing: rotation into the CLsat1 anchor frame (Methods 3.4) and the frame check."""
import os
import random

import pytest

from satpipe import anchor as A

CONS = os.path.join(os.path.dirname(__file__), "..", "satpipe", "stages", "02_seed_expansion_scan", "consensus")
MOTIF = "AAAGAAG"


def _cons(fam):
    return "".join(A.read_fasta(os.path.join(CONS, fam + ".cons.fa")).values())


@pytest.fixture(scope="module")
def frame():
    return A.anchor_frame(_cons("CLsat1"), MOTIF)


def test_shipped_consensuses_are_in_frame(frame):
    assert frame.startswith(MOTIF)
    for fam in ("CLsat1", "CLsat2", "CLsat3", "CLsat4"):
        assert A.verify_records({fam: _cons(fam)}, frame, MOTIF) == []


@pytest.mark.parametrize("fam", ["CLsat1", "CLsat2", "CLsat3", "CLsat4"])
@pytest.mark.parametrize("shift", [1, 36, 70, 111])
def test_any_rotation_is_restored(frame, fam, shift):
    s = _cons(fam)
    k = shift % len(s)
    rotated = s[k:] + s[:k]
    assert A.verify_records({"x": rotated}, frame, MOTIF)          # out of frame is detected
    out, rows = A.anchor_records({"x": rotated}, frame, MOTIF)
    assert out["x"] == s
    assert A.verify_records(out, frame, MOTIF) == []


def test_gaps_ignored_in_check(frame):
    s = _cons("CLsat2")
    gapped = s[:10] + "--" + s[10:]
    assert A.verify_records({"g": gapped}, frame, MOTIF) == []


def test_unplaceable_sequence_fails(frame):
    random.seed(1)
    junk = "".join(random.choice("ACGT") for _ in range(150))
    with pytest.raises(A.AnchorError):
        A.anchor_records({"junk": junk}, frame, MOTIF)


def test_motif_absent_from_anchor_is_an_error():
    with pytest.raises(A.AnchorError):
        A.anchor_frame(_cons("CLsat1"), "CCCCCCCCCC")


def test_no_motif_uses_anchor_as_given():
    a = _cons("CLsat1")
    f = A.anchor_frame(a)
    assert f == a
    s = _cons("CLsat3")
    out, _ = A.anchor_records({"x": s[50:] + s[:50]}, f)
    assert A.verify_records(out, f) == []
