#!/usr/bin/env bash
set -euo pipefail

# 1. Path Setup & Versioning
STAMP="${RUN_DATE:-$(date +%Y%m%d)}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
BASE_DIR="$PROJECT_ROOT/usuv_data"

# Results are now inside a dated subfolder
RESULTS_DIR="$BASE_DIR/results/$STAMP"
# We now use the Master Metadata that was updated in the fetch step
METADATA="$BASE_DIR/curated/usuv_master_metadata.csv"

log() { echo "$(date +"[%Y-%m-%d %H:%M:%S]") [TIMETREE] $*" | tee -a "$BASE_DIR/usuv_pipeline_$STAMP.log"; }

log "=== Starting TreeTime Analysis for Version $STAMP ==="

# 2. Locate Versioned Files
# These names match exactly what analysis.sh produces
ALIGNED_FASTA="$RESULTS_DIR/USUV_Aligned_$STAMP.fasta"
ML_TREE="$RESULTS_DIR/USUV_IQTree_$STAMP.treefile"
TIMETREE_OUT="$RESULTS_DIR/timetree"

# 3. Pre-flight Checks
if [ ! -f "$ML_TREE" ]; then
    log "ERROR: ML tree missing at $ML_TREE. Run --tree first."
    exit 1
fi

if [ ! -f "$ALIGNED_FASTA" ]; then
    log "ERROR: Aligned fasta missing at $ALIGNED_FASTA."
    exit 1
fi

# 4. Clean Metadata (Handle M1/Windows carriage return issues)
log "Sanitizing metadata carriage returns..."
perl -i -pe 's/\r//g' "$METADATA"

# 5. Run TreeTime

log "Estimating molecular clock with TreeTime..."
treetime --aln "$ALIGNED_FASTA" \
         --tree "$ML_TREE" \
         --dates "$METADATA" \
         --name-column "taxa_name" \
         --date-column "collection_date" \
         --outdir "$TIMETREE_OUT" \
         --confidence \
         --reroot "best"

# 6. Final Status Check
if [ -d "$TIMETREE_OUT" ]; then
    log "Success: Time-resolved tree generated in $TIMETREE_OUT"
    # Detailed terminal output you like:
    log "Molecular clock results saved to: $TIMETREE_OUT/molecular_clock.txt"
else
    log "Error: TreeTime failed to create output directory."
    exit 1
fi

log "=== TreeTime Complete for $STAMP ==="