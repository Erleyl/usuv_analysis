#!/usr/bin/env bash
set -euo pipefail

STAMP="${RUN_DATE:-$(date +%Y%m%d)}"
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/usuv_data"
INPUT="$BASE_DIR/results/$STAMP/USUV_QC_Passed_$STAMP.fasta"
OUTPUT="$BASE_DIR/results/$STAMP/USUV_Aligned_$STAMP.fasta"

echo "[ALIGN] Running MAFFT..."
mafft --auto --thread 6 --reorder "$INPUT" > "$OUTPUT"