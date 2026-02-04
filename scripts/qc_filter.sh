#!/usr/bin/env bash
set -euo pipefail

STAMP="${RUN_DATE:-$(date +%Y%m%d)}"
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/usuv_data"
INPUT="$BASE_DIR/curated/USUV_Merged_$STAMP.fasta"
OUTPUT="$BASE_DIR/results/$STAMP/USUV_QC_Passed_$STAMP.fasta"
mkdir -p "$(dirname "$OUTPUT")"

echo "[QC] Filtering: Min Length ${MIN_LEN}bp, Max Ns ${MAX_NS}%"
seqkit seq -m "$MIN_LEN" "$INPUT" | \
seqkit fx2tab -n -g -l | awk -v ns="$MAX_NS" '$4 <= ns' | \
cut -f 1 | seqkit grep -f - <(seqkit seq -m "$MIN_LEN" "$INPUT") > "$OUTPUT"

echo "[QC] Done. Sequences remaining: $(grep -c ">" "$OUTPUT")"