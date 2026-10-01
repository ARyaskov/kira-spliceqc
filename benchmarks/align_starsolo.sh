#!/usr/bin/env bash
# STARsolo alignment template producing all three input levels.
# Usage: align_starsolo.sh <fastq_dir> <star_index> <whitelist.txt> <out_dir> [threads]
# Expects 10x-style R1 (barcode+UMI) / R2 (cDNA) FASTQs: <fastq_dir>/*_R1*.fastq.gz, *_R2*.fastq.gz
set -euo pipefail
FQ=${1:?fastq dir}; IDX=${2:?STAR index}; WL=${3:?whitelist}; OUT=${4:?out dir}; T=${5:-16}
R1=$(ls "$FQ"/*_R1*.fastq.gz | paste -sd, -)
R2=$(ls "$FQ"/*_R2*.fastq.gz | paste -sd, -)
mkdir -p "$OUT"
STAR --runThreadN "$T" --genomeDir "$IDX" \
  --readFilesIn "$R2" "$R1" --readFilesCommand zcat \
  --soloType CB_UMI_Simple --soloCBwhitelist "$WL" \
  --soloUMIlen 12 --soloBarcodeReadLength 0 \
  --soloFeatures Gene Velocyto SJ \
  --soloCellFilter EmptyDrops_CR \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix "$OUT/"
echo "inputs: $OUT/Solo.out/Gene/filtered (L0), $OUT/Solo.out/Velocyto/filtered (L1), $OUT/Solo.out/SJ/raw (L2)"
