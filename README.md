# CLsat workflow → satpipe

A reusable satellite-DNA locus pipeline: one config file, one driver, stage-level inputs and
outputs, and the same viewer ([Toki-bio/CLsat-viewer](https://github.com/Toki-bio/CLsat-viewer))
served from a run's output directory. It is the CLsat/DarSat analysis code from
*"Locus-resolved comparative analysis of minisatellite-scale satellite DNA in bisexual and hybrid
parthenogenetic Darevskia lizards"*, being made runnable on other genomes and other satellite
families. The published Darevskia catalogue is the regression test: the Darevskia config must
reproduce it to the tolerances in [PIPELINE.md](PIPELINE.md).

**Manuscript code deposit.** The code exactly as deposited for the manuscript is tag
[`manuscript-deposit`](../../tree/manuscript-deposit) (the former `main`). Its reproduction
report and every later change that can alter published numbers are in [PIPELINE.md](PIPELINE.md),
which remains the scientific record. The same scripts now live under `satpipe/stages/` with their
original directory names; `git log --follow` gives their history.

`satpipe` is a working name. The final name and the licence are owner decisions
(see [HANDOFF.md](HANDOFF.md), "open decisions").

## Status (work-plan step 3 of 10, 2026-10-15)

| part | state |
|---|---|
| stage scripts | moved with `git mv` from the top level into `satpipe/stages/`, contents unchanged ([PROVENANCE.md](PROVENANCE.md)); superseded code in `legacy/`; the scripts still carry `/staging/...` paths (step 4 strips them) |
| `config.yaml` schema | done: `satpipe/schema.py`, reference in [docs/config.md](docs/config.md); every hard-coded parameter is a key with the published value as default |
| `satpipe check` | done: schema, cross-key checks, input files (genomes, TRF `.dat`, consensus length vs `mono`, curated cluster table, sub-subfamily consensuses, manual loci), external tools with versions |
| `satpipe run`, `satpipe viewer` | not implemented (steps 5-8); they exit 2 with a message |
| Darevskia config | `configs/darevskia_clsat.yaml`; schema-valid; **file paths not yet verified on the server** (see the comment block at its top) |
| tests | `pytest`: 19 tests on a synthetic project, no external tools needed |

## Install

```bash
git clone https://github.com/Toki-bio/CLsat_workflow && cd CLsat_workflow
python3 -m pip install -e ".[test]"      # Python >= 3.10, PyYAML
satpipe --version
python3 -m pytest
```

The pipeline itself needs TRF, VSEARCH, MAFFT, FASTA (`ssearch36`), samtools, bedtools, GNU
parallel, `gawk` (not mawk) and, for optional stages, EMBOSS and BLAST+. `satpipe check` reports
which are missing and which stages they block.

## Use

```bash
satpipe check configs/darevskia_clsat.yaml          # schema + files + tools; exit 1 on any failure
satpipe check my.yaml --no-tools                    # author a config on a machine without the tools
satpipe check my.yaml --no-files --json             # machine-readable
satpipe config my.yaml                              # normalised config: defaults filled, paths resolved
satpipe config --schema-doc > docs/config.md        # regenerate the reference
satpipe stages                                      # stage contracts and wiring status
```

Minimal config ([configs/example_minimal.yaml](configs/example_minimal.yaml)):

```yaml
run_name: my_satellite_run
out_dir: runs/my_satellite_run
genomes:
  spA: /data/genomes/species_A.fna
  spB: /data/genomes/species_B.fna
families:
  SatX: {consensus: ref/SatX.cons.fa, mono: 180}
family_assignment: ref/cluster_family.tsv   # curated cluster -> family table, filled after stage 08
```

Everything not set takes the published Darevskia value, so the behaviour of a minimal config is
the published method applied to new data. Relative paths resolve against the config file's
directory. Unknown keys are errors.

## Layout

```
satpipe/            Python package: cli, config, schema, check, tools, stages registry
satpipe/stages/     the stage scripts, by original directory (01_ ... 07_); consensuses and
                    the curated cluster_family.tsv are in 02_seed_expansion_scan/
configs/            darevskia_clsat.yaml (regression target), example_minimal.yaml
ref/darevskia/      manual_loci.bed (L541, the hand-added locus)
legacy/             superseded code kept for provenance (not current)
viewer/             pointer to the viewer repository; data-product -> producer table
docs/config.md      generated schema reference
docs/files.md       intermediate files per stage (columns, coordinate frames)
tests/              pytest suite (synthetic data)
PIPELINE.md         scientific record: reproduction report, changes after it
HANDOFF.md          mission, conventions, target design, ordered work plan, do-not list
PROVENANCE.md       what moved where; the manuscript-deposit tag
```

## Conventions that must not change

From HANDOFF.md, in short: BED and JSON loci are 0-based half-open; TRF, `samtools faidx` and
ssearch36 `-m 8C` coordinates are 1-based; Δ = assigned-family bits − strongest other family;
`LC_ALL=C` everywhere a sort decides locus numbering; `gawk` required; never change a default
that alters the published catalogue (`short_segments: drop`, plus-strand canonicalisation,
coarse-probe test, scaffold-end handling) — add an option with the published value as default.
