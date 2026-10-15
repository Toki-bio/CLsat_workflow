# Claude Code: read HANDOFF.md first

This repository is the reusable-tool refactor of the CLsat/DarSat satellite-DNA workflow
(work-plan step 3 done on 2026-10-15; steps 4-10 remain). Before editing anything:

1. Read `HANDOFF.md` in full: mission, current state, stage contracts, conventions that must
   not change, target design, ordered work plan, environment gotchas, do-not list.
2. Read `PROVENANCE.md`: what moved where, the `manuscript-deposit` tag, what is in `legacy/`.
   `PIPELINE.md` stays the scientific record; any change that alters published numbers goes there.
3. Run `python3 -m pytest` and `satpipe check configs/darevskia_clsat.yaml --no-files --no-tools`.

Rules in short:

- Stage scripts under `satpipe/stages/` are the tested units of work. Wrap them; do not rewrite
  them in another language until the driver reproduces the Darevskia catalogue with them.
- Never change a default that alters the published numbers. Add an option; keep the published
  value as its default; the schema in `satpipe/schema.py` is the single place parameters live.
- Coordinates: BED/JSON 0-based half-open; TRF/faidx/ssearch 1-based. `LC_ALL=C`. `gawk`.
- Anything that needs genomes or the aligners runs on the analysis server (`/staging`, conda
  env `bioinfo`), not here. Test locally with synthetic data (`tests/`).
- After `satpipe/schema.py` changes, regenerate `docs/config.md` with `satpipe config --schema-doc`.
- Record any change that alters published numbers in `CLsat_workflow/PIPELINE.md` as well.

Steps 3 and 4 are done (see HANDOFF.md progress log). Next is work-plan step 5: the driver
(`satpipe run`, plain Python stage runner, owner's decision) wiring stages 03-14 with DONE
markers and a manifest; check: the Darevskia run from the `.dat` files reproduces the
`cross_species_centroids.fa` md5 and 540 MASTER loci. The toy equivalence harness
(`tests/toy/run_toy_equivalence.sh`, needs the tools; therioserver env `sinederella` has them)
is the regression test for any further script change.
