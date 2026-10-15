#!/usr/bin/env bash
# Re-arbitration of a corrected catalogue = Step 14 run again on the corrected MASTER BED with
# the rescan bits of every family. In the published run this file was a copy of
# arbitration_patched.sh with ten path lines patched by sed (run_arbitration_540.sh, now in
# legacy/published_run/); since the script takes its inputs as arguments there is nothing to
# patch, so this is a thin wrapper kept under its historical name.
#
# Usage: arbitration_540.sh <arbitration_inputs.tsv> <MASTER.merged.bed> <out_dir>
#   (same arguments and environment as 03_master_arbitration/arbitration_patched.sh)
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec bash "$HERE/../../03_master_arbitration/arbitration_patched.sh" "$@"
