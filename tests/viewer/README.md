# Viewer tests

- `compare_embedded.py` and `make_config_from_viewer_dir.py` (list rebuild check and the one-time migration of the
  live Darevskia site) used the Python driver modules and are at tag `python-driver-final`.
- `cdp_test.js <dir> <port>`: serves `<dir>` with `python -m http.server`, opens the page in headless
  Chrome (Windows path to Chrome inside), waits for the data load, exercises the species buttons, a locus
  click, the three inspection lists and the sub-subfamily filter, reports console errors and writes
  `screenshot.png`. Run after any change to `satpipe/viewer_html/index.html`.
- `darevskia_viewer_meta.json`: the constants that used to be embedded in the Darevskia page.
