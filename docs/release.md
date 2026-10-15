# Releasing

1. Decide the two things still open in HANDOFF.md: the tool name (rename this repository or keep
   `CLsat_workflow`; GitHub redirects the old URL) and the licence (write `LICENSE`, set `license`
   in `CITATION.cff` and `pyproject.toml`).
2. Bump `version` in `pyproject.toml`, `satpipe/__init__.py` and `CITATION.cff`; set `date-released`.
3. Make sure CI is green (`.github/workflows/ci.yml`: unit tests without tools; toy equivalence and
   driver tests in the conda environment).
4. Tag: `git tag -a v0.1.0 -m "satpipe v0.1.0" && git push origin v0.1.0`, then create the GitHub release
   from the tag.
5. Zenodo: enable the repository under https://zenodo.org/account/settings/github/ before the release
   is published; Zenodo archives the release and mints the DOI. Put the DOI badge in README.md and the
   `doi` field in `CITATION.cff` afterwards.
6. The manuscript's data-availability statement should cite the tag `manuscript-deposit` (the code as
   deposited) and the release DOI (the tool).

Regenerating the published viewer after a data change in Toki-bio/CLsat-viewer:

    python tests/viewer/make_config_from_viewer_dir.py <viewer clone> <viewer clone>/viewer_meta.json

and commit `viewer_config.json`. After a change to the page, re-vendor it:

    python tools/patch_viewer_config.py <pristine v2.25 index.html> satpipe/viewer_html/index.html
    node tests/viewer/cdp_test.js <viewer dir> 18931
