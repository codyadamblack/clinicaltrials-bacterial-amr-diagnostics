#!/usr/bin/env bash
set -Eeuo pipefail

if [[ "$#" -lt 2 ]]; then
    echo "Usage: $0 INPUT_DIR OUTPUT_DIR [additional screening arguments...]" >&2
    exit 2
fi

INPUT_DIR="$1"
OUTPUT_DIR="$2"
shift 2

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

python "$ROOT/screening/screen_bacterial_amr_diagnostics_v3_2_5.py" \
    --input-dir "$INPUT_DIR" \
    --output-dir "$OUTPUT_DIR" \
    "$@"
