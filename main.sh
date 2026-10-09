#!/usr/bin/env bash
set -euo pipefail

# 1. GLOBAL VERSIONING & DEFAULTS
export RUN_DATE="${RUN_DATE:-$(date +%Y%m%d)}"
export FORCE=0
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PROJECT_ROOT="$SCRIPT_DIR"

# QC Thresholds (Defaults can be overridden via command line)
# Defaults follow Emmanuelle Munger's rule (240729_PrepSeqData.R): longer than
# 8800 nt and fewer than 2200 N (about 20% N for an 11 kb genome).
# Before: MIN_LEN=10000, MAX_NS=10.0. Keep in sync with scripts/qc_filter.sh.
export MIN_LEN=${MIN_LEN:-8800}
export MAX_NS=${MAX_NS:-20.0}

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
    echo "  Example: MIN_LEN=5000 MAX_NS=5.0 ./main.sh --qc"
    echo ""
    echo "OPTIONS:"
    echo "  --dry-run     SYSTEM CHECK: Verifies tools and project structure."
    echo "  --fetch       FETCH: Incremental NCBI update (Master Database)."
    echo "  --lab         LAB: Prepare lab + partner (e.g. RIVM) sequences for the merge (R)."
    echo "                Settings: LAB_SHEET, LAB_EXTERNAL, LAB_R_ARGS (see scripts/prepare_lab.sh)."
    echo "  --merge       MERGE: Master + Lab data & Exclude Blacklist. Also writes the tip-label table."
    echo "  --qc          QC: Filter by Length ($MIN_LEN) and Ns ($MAX_NS%)."
    echo "  --align       ALIGN: MAFFT Multiple Sequence Alignment."
    echo "  --tree        TREE: IQ-TREE Maximum Likelihood Phylogeny."
    echo "  --local-tree  TREE (LOCAL ONLY): IQ-TREE on local-lab-curated samples only."
    echo "  --time        TIME: TreeTime Molecular Clock Analysis."
    echo "  --local-time  TIME (LOCAL ONLY): TreeTime on the local-lab-curated tree."
    echo "  --labels      LABELS: Tip-label table (TSV) for PearTree. Made automatically by --merge;"
    echo "                use this to remake it (e.g. with another LABEL_TEMPLATE)."
    echo "                Settings: LABEL_TEMPLATE, LABEL_LOWER (see scripts/tip_labels.sh)."
    echo "  --force       FORCE: redo steps already run (see below). Combine with other options."
    echo "  --all         FULL RUN: Executes all steps sequentially."
    echo ""
    echo "EXAMPLES:"
    echo "  bash main.sh --fetch --lab --merge --qc --align --tree --time"
    echo "  bash main.sh --lab --merge --qc --align --local-tree  # lab data changed, no new NCBI fetch"
    echo "  LAB_EXTERNAL=none bash main.sh --lab                  # without partner (RIVM) data"
    echo "  bash main.sh --labels                                 # tip-label table for the PearTree view of this run"
    echo "  bash main.sh --force --merge --qc --align --tree      # redo a finished run"
    echo "  RUN_DATE=20260730 bash main.sh --force --tree         # redo an older run"
    echo "  FETCH_SINCE=2026-06-01 bash main.sh --force --fetch   # re-download since a date"
    echo "========================================================="
    exit 0
}

[ $# -eq 0 ] && usage

START_TIME=$(date +%s)

# Pre-scan for --force so it applies to every step regardless of its position
for arg in "$@"; do
    [[ "$arg" == "--force" ]] && FORCE=1
done
if [[ "$FORCE" == "1" ]]; then
    log "FORCE mode: --fetch re-downloads from ${FETCH_SINCE:-1900-01-01}; IQ-TREE runs with -redo."
    export FETCH_SINCE="${FETCH_SINCE:-1900-01-01}"
fi

for arg in "$@"; do
    case $arg in
        --force)
            : # handled in the pre-scan above
            ;;
        --dry-run)
            log "DRY RUN: Verifying software dependencies..."
            for cmd in datasets mafft iqtree treetime seqkit python3 Rscript; do
                command -v $cmd &> /dev/null && echo "  [OK] $cmd" || echo "  [!!] $cmd MISSING"
            done
            ;;
        --fetch)
            log "STEP 1: Starting Incremental Update from NCBI..."
            bash "$PROJECT_ROOT/scripts/update_usuv.sh"
            ;;
        --lab)
            log "STEP 1b: Preparing lab + partner sequences (R: clean, drop what is already on NCBI)..."
            bash "$PROJECT_ROOT/scripts/prepare_lab.sh"
            ;;
        --merge)
            log "STEP 2: Merging Lab Data and Applying Exclusion List..."
            bash "$PROJECT_ROOT/scripts/merge_lab.sh"
            # the tip-label table (for PearTree) only needs the merged metadata: make it right away
            bash "$PROJECT_ROOT/scripts/tip_labels.sh" || log "WARNING: tip-label table could not be made (pipeline continues)."
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
        --local-tree)
            log "STEP 5b: Running IQ-TREE on local-lab-curated samples only..."
            bash "$PROJECT_ROOT/scripts/local_lab_tree.sh"
            ;;
        --time)
            log "STEP 6: Running TreeTime Molecular Clock..."
            bash "$PROJECT_ROOT/scripts/timetree.sh"
            ;;
        --local-time)
            log "STEP 6b: Running TreeTime on local-lab-curated tree..."
            bash "$PROJECT_ROOT/scripts/local_lab_timetree.sh"
            ;;
        --labels)
            log "STEP 7: Making the tip-label table (for PearTree)..."
            bash "$PROJECT_ROOT/scripts/tip_labels.sh"
            ;;
        --all)
            log "FULL PIPELINE INITIATED"
            bash "$PROJECT_ROOT/scripts/update_usuv.sh"
            bash "$PROJECT_ROOT/scripts/prepare_lab.sh"
            bash "$PROJECT_ROOT/scripts/merge_lab.sh"
            bash "$PROJECT_ROOT/scripts/tip_labels.sh" || log "WARNING: tip-label table could not be made (pipeline continues)."
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
