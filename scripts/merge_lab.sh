#!/usr/bin/env bash
set -euo pipefail

# 1. Environment & Path Setup
STAMP="${RUN_DATE:-$(date +%Y%m%d)}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

BASE_DIR="$PROJECT_ROOT/usuv_data"
CURATED_DIR="$BASE_DIR/curated"
EXCLUDE_FILE="$PROJECT_ROOT/exclude_accessions.txt"

# Master file from fetch step and the new dated snapshot
MASTER_FASTA="$CURATED_DIR/usuv_master.fasta"
FINAL_FASTA="$CURATED_DIR/USUV_Merged_$STAMP.fasta"
LOCAL_LAB_FASTA="$BASE_DIR/local_lab/internal_genomes.fasta"

mkdir -p "$BASE_DIR/local_lab"
log() { echo "$(date +"[%Y-%m-%d %H:%M:%S]") [MERGE] $*" | tee -a "$BASE_DIR/usuv_pipeline_$STAMP.log"; }

# 2. Start building the snapshot
log "Combining Master database with local lab data for version $STAMP..."

if [ ! -f "$MASTER_FASTA" ]; then
    log "ERROR: Master database not found at $MASTER_FASTA. Run --fetch first."
    exit 1
fi

# Copy master to the new dated file
cat "$MASTER_FASTA" > "$FINAL_FASTA"

# 3. Add Lab Data
if [ -f "$LOCAL_LAB_FASTA" ]; then
    cat "$LOCAL_LAB_FASTA" >> "$FINAL_FASTA"
    log "Lab data merged successfully."
else
    log "No lab data found at $LOCAL_LAB_FASTA. Using NCBI only."
fi

# 4. Apply Exclusion List (The "Bouncer" Step)
if [ -f "$EXCLUDE_FILE" ]; then
    log "Checking exclusion list..."
    # count IDs in file
    COUNT=$(grep -c '^' "$EXCLUDE_FILE" || echo 0)
    if [ "$COUNT" -gt 0 ]; then
        log "Removing $COUNT blacklisted accessions from this run..."
        seqkit grep -v -f "$EXCLUDE_FILE" "$FINAL_FASTA" -o "$FINAL_FASTA.tmp"
        mv "$FINAL_FASTA.tmp" "$FINAL_FASTA"
    fi
else
    log "No exclusion list found at $EXCLUDE_FILE. Skipping blacklist step."
fi

log "Final snapshot created: $FINAL_FASTA"