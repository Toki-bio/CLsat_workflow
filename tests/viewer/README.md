# Viewer tests

- `compare_embedded.py <viewer_dir>`: rebuilds ORTHO_LOCI / FLAGGED_LOCI / FULL_TABLE from the data
  files with `satpipe.viewer_lists` and compares them with the lists embedded in a pre-v3 `index.html`
  (the v2.25 audit build). ORTHO_LOCI and FLAGGED_LOCI rebuild identically; FULL_TABLE differs only in
  row order and in the `both` flag (see PROVENANCE.md, "both-sided").
- `make_config_from_viewer_dir.py <viewer_dir> darevskia_viewer_meta.json`: writes `viewer_config.json`
  for an already-deployed viewer directory (migration of the live Darevskia site to the config-driven page).
- `cdp_test.js <dir> <port>`: serves `<dir>` with `python -m http.server`, opens the page in headless
  Chrome (Windows path to Chrome inside), waits for the data load, exercises the species buttons, a locus
  click, the three inspection lists and the sub-subfamily filter, reports console errors and writes
  `screenshot.png`. Run after any change to `satpipe/viewer_html/index.html`.
- `darevskia_viewer_meta.json`: the constants that used to be embedded in the Darevskia page.
