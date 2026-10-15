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

## Status (work-plan steps 1-8 done, 9-10 in progress, 2026-10-15)

| part | state |
|---|---|
| stage scripts | in `satpipe/stages/`, parameterised: inputs from arguments and small tables (`genomes.tsv`, `families.tsv`, `members.tsv`), parameters from environment variables with the published values as defaults; no server path, species code, family name or monomer length in code ([PROVENANCE.md](PROVENANCE.md)); equivalence with the original scripts verified on a synthetic project (`tests/toy/`) |
| `config.yaml` schema | done: `satpipe/schema.py`, reference in [docs/config.md](docs/config.md); every hard-coded parameter is a key with the published value as default |
| `satpipe check` | done: schema, cross-key checks, input files, external tools with versions |
| `satpipe run` | done: 20 stages from TRF to flank orthology, resumable (DONE markers), `manifest.json`; `--from/--to/--only/--force` |
| `satpipe viewer` | done: assembles `<run>/viewer/` (data files + `viewer_config.json` + the config-driven page) |
| Darevskia regression | running on therioserver from the six GenBank assemblies (`configs/darevskia_clsat_therioserver.yaml`); see HANDOFF.md progress log |
| second test case | Amazon molly (`tests/molly/run_molly.sh`), running |
| tests | `pytest` (no tools); `tests/toy/` equivalence + driver tests (tools); `tests/viewer/` list rebuild + headless-Chrome page test; CI in `.github/workflows/ci.yml` |

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
satpipe run my.yaml                                 # all stages; re-run = resume (DONE markers)
satpipe run my.yaml --to 08_cross                   # stop after cross-clustering, inspect cross/cluster_report.tsv,
                                                    #   write the curated cluster -> family table, then run again
satpipe run my.yaml --from 13_rescan --force        # redo from a stage
satpipe viewer runs/my_run                          # -> runs/my_run/viewer/ : serve it with any static server
satpipe config my.yaml                              # normalised config: defaults filled, paths resolved
satpipe config --schema-doc > docs/config.md        # regenerate the reference
satpipe run --list-stages my.yaml                   # stage ids
```

A run from scratch: `satpipe run` runs TRF on every genome (stage 02; give `inputs.trf_dat` to skip it),
extracts and clusters motifs, cross-clusters them with the family references prepended, and stops at
stage 09 if the curated `family_assignment` table does not match the clusters it found. Look at
`<run>/cross/cluster_report.tsv` (every cluster with its size, centroid and best family hit), write
the table, run again. Everything after that is automatic: seeds, expansion, MASTER, rescan,
arbitration, catalogue (+ manual loci), regions, terminal monomers, sub-subfamilies (if consensuses
are given), flank orthology, and `satpipe viewer` for the page.

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
