# lib/satpipe - the satpipe driver (bash + gawk)

`bin/satpipe` sources these files. One file per stage: `st_<id>.sh` defines `plan_<id>` (inputs + params for the DONE
fingerprint) and `st_<id>` (the commands). `schema.tsv` is the parameter schema (every config key, type, default, doc);
edit it directly and regenerate `docs/config.md` with `satpipe config --schema-doc`.

Comments of the form `run.py:N`, `cli.py:N`, `viewer.py:N` point at the Python driver this code replaced, at tag
`python-driver-final` (the port copied its behaviour line by line; the equivalence harness is there too, `tests/port/`).
