# 07 — post-catalogue checks and figure regeneration

Everything here runs **after** the main pipeline (steps 01–06) has produced the
626-locus catalogue. None of it modifies that catalogue, with one documented
exception: `singleton_screen.sh` is what recovered **L541**, the 627th segmented
locus, which the pipeline cannot produce (see the note in `../PIPELINE.md`).

These scripts back specific statements in the manuscript. They are deposited so
those statements can be checked, not because the main pipeline depends on them.

| script | what it establishes | where it is used |
|---|---|---|
| `singleton_screen.sh` | screens all 979 singleton loci against the four CLsat consensuses | Methods §3.1; recovered L541; identified the seven uncatalogued arrays |
| `candidate_scan.py` | four-way per-window scan + arbitration for those seven | Discussion (small, diverged arrays); viewer flag `CANDIDATE_UNCATALOGUED` |
| `candidate_subsubfamily.py` | sub-subfamily call for each candidate | as above |
| `candidate_flank5k.py` | rescans each candidate with 5 kb of context | shows the arrays are genuinely short, not seed-truncated |
| `imotaq_compare.sh` | compares all GenBank IMO-TaqI records against the 39 DarSat consensuses | Results/Discussion: DarSat2 **is** IMO-TaqI |
| `make_fig2_from_data.py` | renders Figure 2 from the deposited catalogue | Figure 2 |

## Why a singleton screen exists at all

Per-assembly clustering (step 06) runs `vsearch --cluster_fast --id 0.90` at its
default `--strand plus`. A monomer and its reverse complement therefore cluster
separately. Across the six assemblies the 7,546 candidate seed loci form 1,593
clusters: 614 holding two or more loci (6,567 loci in total) and **979 singleton
loci**. Only the aligned members of the multi-locus clusters enter the
cross-assembly pass, so a singleton locus never becomes a seed — regardless of how
strong its CLsat signal is.

The screen asks what that cost. Of the 979, **14** score above the 50-bit noise
floor, all of them above 77 bits:

* **7** already lie inside catalogued loci, swallowed by a neighbouring seed's
  boundary expansion;
* **1** is L541 (mix, `JAWWNH010000226.1:1–43,510`, ~298 monomers, CLsat3 at 127.8
  bits) — added to the catalogue by hand;
* **7** remain outside it. They are short (1.5–2.0 kb) and diverged (66.7–80.8 %
  identity to their nearest sub-subfamily consensus, against a catalogue median of
  97.4 %). They are **excluded from every number in the manuscript** and appear only
  in the viewer, flagged `CANDIDATE_UNCATALOGUED` with ids `C1`–`C7`.

## Note on Figure 2

Figure 2 was originally rendered from a hardcoded count matrix, which silently went
stale when L541 was added (`CLsat3.1/mix 2→3`) and a val locus was reassigned
(`CLsat1.5/val 1→0`) — the figure showed mix = 128, val = 80 against the catalogue's
129 and 79. Both versions summed to 627, so total-based checks never caught it.
`make_fig2_from_data.py` derives the matrix from `MASTER.loci.segmented.bed` and
`locus_subfamily.tsv` and cross-checks per-assembly totals against an independent
count of the BED, so it cannot drift again.

## Paths

These were run on the analysis server and carry absolute paths under
`/staging/Sergei/Lacertidae/Darevskia/`. They are recorded as run, not rewritten as
portable tools; adjust the paths at the top of each script to re-execute.
`ssearch36` is from the FASTA package, `blastn`/`makeblastdb` from BLAST+.
