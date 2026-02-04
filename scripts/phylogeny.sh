#!/usr/bin/env bash
set -euo pipefail

STAMP="${RUN_DATE:-$(date +%Y%m%d)}"
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/usuv_data"
INPUT="$BASE_DIR/results/$STAMP/USUV_Aligned_$STAMP.fasta"
PRE="$BASE_DIR/results/$STAMP/USUV_IQTree_$STAMP"

echo "[TREE] Running IQ-TREE (MFP + 1000 Bootstraps)..."
iqtree -s "$INPUT" -m MFP -bb 1000 -alrt 1000 -nt AUTO -pre "$PRE"