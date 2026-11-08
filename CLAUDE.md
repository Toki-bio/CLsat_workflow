# Claude Code: read HANDOFF.md first

This repository is the reusable-tool refactor of the CLsat/DarSat satellite-DNA workflow
(work-plan step 3 done on 2026-10-15; steps 4-10 remain). Before editing anything:

1. Read `HANDOFF.md` in full: mission, current state, stage contracts, conventions that must
   not change, target design, ordered work plan, environment gotchas, do-not list.
2. Read `PROVENANCE.md`: what moved where, the `manuscript-deposit` tag, what is in `legacy/`.
   `PIPELINE.md` stays the scientific record; any change that alters published numbers goes there.
3. Run `bash tests/run_all.sh` and `bin/satpipe check configs/darevskia_clsat.yaml --no-files --no-tools` (POSIX host).

Rules in short:

- Stage scripts under `satpipe/stages/` are the tested units of work. Wrap them; do not rewrite
  them in another language until the driver reproduces the Darevskia catalogue with them.
- Never change a default that alters the published numbers. Add an option; keep the published
  value as its default; the schema in `lib/satpipe/schema.tsv` is the single place parameters live.
- Coordinates: BED/JSON 0-based half-open; TRF/faidx/ssearch 1-based. `LC_ALL=C`. `gawk`.
- Anything that needs genomes or the aligners runs on the analysis server (`/staging`, conda
  env `bioinfo`), not here. Test locally with synthetic data (`tests/`).
- After `lib/satpipe/schema.tsv` changes, regenerate `docs/config.md` with `bin/satpipe config --schema-doc`.
- The driver is bash + gawk (`bin/satpipe`, `lib/satpipe/`) since 2026-11-09; the Python driver is at tag
  `python-driver-final`. Do not reintroduce a Python driver.
- Record any change that alters published numbers in `CLsat_workflow/PIPELINE.md` as well.

Work-plan state: HANDOFF.md progress log. The toy equivalence harness (`tests/toy/run_toy_equivalence.sh`,
needs the tools; therioserver env `sinederella` has them) is the regression test for any stage-script change.
