#!/usr/bin/env bash
set -Eeuo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="$(mktemp)"
trap 'rm -f "$OUT"' EXIT

python "$ROOT/screening/screen_bacterial_amr_diagnostics_v3_2_5.py" \
    --self-test \
    | tee "$OUT"

grep -Fq "Version 3.2.5 self-tests: PASS" "$OUT"

echo "SCREENING SELF-TEST WRAPPER: PASS"
