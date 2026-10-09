#!/usr/bin/env bash
# ==============================================================================
# Script Name: tip_labels.sh              (main.sh step: --labels)
# Description: Writes a tip-label table for the run, to be loaded into PearTree (or
#              FigTree / ggtree) as annotations:
#                  usuv_data/results/<RUN_DATE>/USUV_tip_labels_<RUN_DATE>.tsv
#              The first column is taxa_name (= the tip names of the trees), followed by
#              host, date, source, project, public, lab_id, country, location, year and one label_<style> column per label
#              style. The trees and FASTA files are NOT changed.
#              Only the merged metadata is needed. main.sh runs this step automatically at the
#              end of --merge (and --all); run --labels to remake it, e.g. with LABEL_TEMPLATE.
#
# Settings (environment variables, all optional):
#   LABEL_TEMPLATE  extra style, e.g. "{source}_{acc}_{country3}_{year}" (column label_custom)
#   LABEL_LOWER=1   lower-case host and country
#   See all styles/fields: python3 scripts/make_tip_labels.py --list-styles
# ==============================================================================
set -euo pipefail

STAMP="${RUN_DATE:-$(date +%Y%m%d)}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
BASE_DIR="$PROJECT_ROOT/usuv_data"
RESULTS_DIR="$BASE_DIR/results/$STAMP"
OUT="$RESULTS_DIR/USUV_tip_labels_$STAMP.tsv"

mkdir -p "$BASE_DIR/logs" "$RESULTS_DIR"
log()  { echo "$(date +"[%Y-%m-%d %H:%M:%S]") [LABELS] $*" | tee -a "$BASE_DIR/usuv_pipeline_$STAMP.log"; }
fail() { log "ERROR: $*"; exit 1; }

command -v python3 >/dev/null 2>&1 || fail "python3 not found."

# Metadata of this run (the one the sequences/trees were built from)
META="$BASE_DIR/curated/$STAMP/USUV_Merged_${STAMP}_metadata.csv"
[[ -f "$META" ]] || META="$BASE_DIR/curated/USUV_Merged_${STAMP}_metadata.csv"   # legacy flat layout
[[ -f "$META" ]] || fail "Merged metadata for $STAMP not found. Run --merge first (or set RUN_DATE)."

ARGS=(--metadata "$META" --out "$OUT")

# Selection audit of the --lab step: tells which NCBI sequences belong to the lab / partners
AUDIT="$BASE_DIR/lab_sequences/USUV_lab_${STAMP}_selection_audit.tsv"
if [[ ! -f "$AUDIT" ]]; then
    AUDIT="$(ls -1 "$BASE_DIR"/lab_sequences/USUV_lab_*_selection_audit.tsv 2>/dev/null | sort | tail -n 1 || true)"
    [[ -n "$AUDIT" ]] && log "Note: no lab audit for $STAMP (run --lab); using $(basename "$AUDIT")."
fi
if [[ -n "$AUDIT" && -f "$AUDIT" ]]; then
    ARGS+=(--lab-audit "$AUDIT")
else
    log "WARNING: no lab selection audit found (run --lab): public sequences of the lab will show as NCBI."
fi
[[ -n "${LABEL_TEMPLATE:-}" ]]   && ARGS+=(--template "$LABEL_TEMPLATE")
[[ "${LABEL_LOWER:-0}" == "1" ]] && ARGS+=(--lower)

python3 "$SCRIPT_DIR/make_tip_labels.py" "${ARGS[@]}" 2>&1 | tee -a "$BASE_DIR/usuv_pipeline_$STAMP.log"

log "Tip-label table for PearTree: $OUT"
log "Load it as annotations next to the tree (first column taxa_name = tip name), then pick a label_* column."
