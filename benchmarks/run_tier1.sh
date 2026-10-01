#!/usr/bin/env bash
# Tier 1: simulate -> run -> validate. Usage: run_tier1.sh <out_dir> [n_cells] [seed]
set -euo pipefail
OUT=${1:?out dir}; N=${2:-2000}; SEED=${3:-24301}
BIN=${KIRA_SPLICEQC:-kira-spliceqc}
mkdir -p "$OUT"
"$BIN" simulate --out "$OUT/sim" --n-cells "$N" --seed "$SEED"
"$BIN" run --input "$OUT/sim" --junctions "$OUT/sim/sj" --out "$OUT/run" --run-mode pipeline
"$BIN" validate --run "$OUT/run" --truth "$OUT/sim/truth.tsv" --out "$OUT/validation.json"
