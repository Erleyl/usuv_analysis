#!/usr/bin/env bash
set -euo pipefail

# 1. GLOBAL VERSIONING & DEFAULTS
export RUN_DATE=$(date +%Y%m%d)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PROJECT_ROOT="$SCRIPT_DIR"

# QC Thresholds (Defaults can be overridden via command line)
export MIN_LEN=${MIN_LEN:-10000}
export MAX_NS=${MAX_NS:-5.0}

# Directory paths
BASE_DIR="$PROJECT_ROOT/usuv_data"
LOG_DIR="$BASE_DIR/logs"
mkdir -p "$LOG_DIR"

log() {
    printf "\033[1;34m$(date +"[%Y-%m-%d %H:%M:%S]") [MAIN] %s\033[0m\n" "$*" | tee -a "$LOG_DIR/pipeline_$RUN_DATE.log"
}

usage() {
    echo -e "\033[1;32m========================================================="
    echo " USUTU VIRUS PHYLOGENY PIPELINE - OPERATING MANUAL"
    echo -e "=========================================================\033[0m"
    echo "Usage: [THRESHOLDS] ./main.sh [OPTION]"
    echo ""
    echo "QC SETTINGS (Current Defaults: $MIN_LEN bp, $MAX_NS % Ns):"
    echo "  Example: MIN_LEN=5000 MAX_NS=2.0 ./main.sh --qc"
    echo ""
    echo "OPTIONS:"
    echo "  --dry-run  SYSTEM CHECK: Verifies tools and project structure."
    echo "  --fetch    FETCH: Incremental NCBI update (Master Database)."
    echo "  --merge    MERGE: Master + Lab data & Exclude Blacklist."
    echo "  --qc       QC: Filter by Length ($MIN_LEN) and Ns ($MAX_NS%)."
    echo "  --align    ALIGN: MAFFT Multiple Sequence Alignment."
    echo "  --tree     TREE: IQ-TREE Maximum Likelihood Phylogeny."
    echo "  --time     TIME: TreeTime Molecular Clock Analysis."
    echo "  --all      FULL RUN: Executes all steps sequentially."
    echo "========================================================="
    exit 0
}

[ $# -eq 0 ] && usage

START_TIME=$(date +%s)

for arg in "$@"; do
    case $arg in
        --dry-run)
            log "DRY RUN: Verifying software dependencies..."
            for cmd in datasets mafft iqtree treetime seqkit python3; do
                command -v $cmd &> /dev/null && echo "  [OK] $cmd" || echo "  [!!] $cmd MISSING"
            done
            ;;
        --fetch)
            log "STEP 1: Starting Incremental Update from NCBI..."
            bash "$PROJECT_ROOT/scripts/update_usuv.sh"
            ;;
        --merge)
            log "STEP 2: Merging Lab Data and Applying Exclusion List..."
            bash "$PROJECT_ROOT/scripts/merge_lab.sh"
            ;;
        --qc)
            log "STEP 3: Running QC Filter (Len: $MIN_LEN, Ns: $MAX_NS%)..."
            bash "$PROJECT_ROOT/scripts/qc_filter.sh"
            ;;
        --align)
            log "STEP 4: Running MAFFT Alignment..."
            bash "$PROJECT_ROOT/scripts/align.sh"
            ;;
        --tree)
            log "STEP 5: Running IQ-TREE Phylogeny..."
            bash "$PROJECT_ROOT/scripts/phylogeny.sh"
            ;;
        --time)
            log "STEP 6: Running TreeTime Molecular Clock..."
            bash "$PROJECT_ROOT/scripts/timetree.sh"
            ;;
        --all)
            log "FULL PIPELINE INITIATED"
            bash "$PROJECT_ROOT/scripts/update_usuv.sh"
            bash "$PROJECT_ROOT/scripts/merge_lab.sh"
            bash "$PROJECT_ROOT/scripts/qc_filter.sh"
            bash "$PROJECT_ROOT/scripts/align.sh"
            bash "$PROJECT_ROOT/scripts/phylogeny.sh"
            bash "$PROJECT_ROOT/scripts/timetree.sh"
            log "FULL PIPELINE COMPLETED"
            ;;
        *)
            usage
            ;;
    esac
done

# Final Runtime Report
END_TIME=$(date +%s)
log "Total execution time: $((END_TIME - START_TIME)) seconds."
