#!/usr/bin/env bash
# Run kira-spliceqc on a STARsolo Gene directory and score it against a truth table.
# Usage: run_dataset.sh <Solo.out/Gene/<subset>> <truth.tsv> <out_dir> [--pair truth:metric:flag ...]
set -euo pipefail
GENE=${1:?Gene dir}; TRUTH=${2:?truth tsv}; OUT=${3:?out dir}; shift 3
BIN=${KIRA_SPLICEQC:-kira-spliceqc}
mkdir -p "$OUT"
# Velocyto/ and SJ/ siblings are auto-detected; metadata.tsv next to the Gene dir gives strata.
"$BIN" run --input "$GENE" --out "$OUT/run" --run-mode pipeline
"$BIN" validate --run "$OUT/run" --truth "$TRUTH" --out "$OUT/validation.json" "$@"
