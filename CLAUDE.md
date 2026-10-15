# Claude Code: read HANDOFF.md first

This repository is the code deposit for the CLsat/DarSat satellite-DNA manuscript (Darevskia lizards)
and the starting point for turning it into a reusable satellite-DNA tool.

- Work on branch `audit-fixes` (or its merge into `main`); `main` before 2026-10-15 has known bugs.
- `HANDOFF.md` holds the mission, current state, stage contracts, conventions that must not change,
  the target design (config + driver + viewer config), the ordered work plan, environment gotchas,
  and a do-not list. Read it in full before editing anything.
- `PIPELINE.md` is the scientific record: reproduction report and "Changes after the reproduction
  report". Any change that alters published numbers must be recorded there.
- Every script still carries `/staging/...` server paths; nothing here runs without the genomes,
  ssearch36, vsearch, MAFFT and TRF. Test with synthetic data (see HANDOFF.md, Verification).
- Companion repository: Toki-bio/CLsat-viewer (branch `audit-fixes`, v2.25).

Conventions in short: BED and JSON loci are 0-based half-open; TRF/faidx/ssearch coordinates are
1-based; Δ = assigned-family bits − strongest other family; keep `LC_ALL=C`; `gawk` required;
never change a default that alters the published catalogue — add an option instead.
