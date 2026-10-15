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

Next step is work-plan step 4 (parameterise the scripts: strip `/staging` constants, the
`declare -A` tables and literal MONO/window widths; feed them from the driver). Check for that
step: `grep -rn '/staging\|JAWWN\|146\b' satpipe/stages` returns nothing that is not a comment.
