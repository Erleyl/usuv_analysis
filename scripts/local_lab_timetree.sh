#!/usr/bin/env bash
#
# local_lab_timetree.sh
# ======================
# Runs TreeTime molecular clock analysis on the local-lab + reference tree
# produced by local_lab_tree.sh (not the full-dataset tree from timetree.sh).
#
# Reads dates from the MERGED metadata (not usuv_master_metadata.csv),
# because local-lab rows only exist there -- the master file is NCBI-only.
#
# Prerequisites: run --local-tree first for the target $RUN_DATE.
#
# Usage:
#   bash local_lab_timetree.sh
#   RUN_DATE=20260730 bash local_lab_timetree.sh
#
set -euo pipefail

STAMP="${RUN_DATE:-$(date +%Y%m%d)}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

BASE_DIR="$PROJECT_ROOT/usuv_data"
CURATED_DIR="$BASE_DIR/curated"
RESULTS_DIR="$BASE_DIR/results/$STAMP"
LOCAL_DIR="$RESULTS_DIR/local_lab_only"

MERGED_META="$CURATED_DIR/$STAMP/USUV_Merged_${STAMP}_metadata.csv"
# Legacy flat layout fallback (runs created before dated folders)
[[ -f "$MERGED_META" ]] || MERGED_META="$CURATED_DIR/USUV_Merged_${STAMP}_metadata.csv"
LOCAL_ALIGNED="$LOCAL_DIR/USUV_LocalLab_WithRefs_Aligned_$STAMP.fasta"
LOCAL_TREE="$LOCAL_DIR/USUV_LocalLab_WithRefs_IQTree_$STAMP.treefile"
TIMETREE_OUT="$LOCAL_DIR/timetree"
META_CLEAN="$LOCAL_DIR/USUV_LocalLab_WithRefs_metadata_clean_$STAMP.csv"

log() { echo "$(date +"[%Y-%m-%d %H:%M:%S]") [LOCAL-TIMETREE] $*" | tee -a "$BASE_DIR/usuv_pipeline_$STAMP.log"; }

log "=== Starting local-lab TreeTime analysis for $STAMP ==="

# 1. Pre-flight checks
if [ ! -f "$LOCAL_TREE" ]; then
    log "ERROR: local-lab ML tree missing at $LOCAL_TREE. Run 'main.sh --local-tree' first."
    exit 1
fi
if [ ! -f "$LOCAL_ALIGNED" ]; then
    log "ERROR: local-lab alignment missing at $LOCAL_ALIGNED. Run 'main.sh --local-tree' first."
    exit 1
fi
if [ ! -f "$MERGED_META" ]; then
    log "ERROR: merged metadata missing at $MERGED_META. Run 'main.sh --merge' first."
    exit 1
fi

# 2. Build a clean, local-scoped copy of the dates file rather than editing
#    the shared merged metadata in place. Also convert literal "unknown"
#    collection_date values to blank, so TreeTime treats them as dates to
#    estimate rather than failing to parse the string "unknown".
log "Preparing dates file (copy + sanitize collection_date)..."
python3 - "$MERGED_META" "$META_CLEAN" <<'PY'
import csv, sys

meta_path, out_path = sys.argv[1], sys.argv[2]

with open(meta_path, newline="", encoding="utf-8") as fh:
    sample = fh.read(4096)
    fh.seek(0)
    delim = ";" if sample.count(";") > sample.count(",") else ","
    rdr = csv.DictReader(fh, delimiter=delim)
    rows = list(rdr)
    fields = rdr.fieldnames or []

n_unknown = 0
for row in rows:
    val = (row.get("collection_date") or "").strip()
    if val.lower() in ("unknown", "nan", "") or val.startswith("["):   # "[2015 TO 2019]" ranges: not readable by TreeTime
        row["collection_date"] = ""
        n_unknown += 1

with open(out_path, "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
    w.writeheader()
    for row in rows:
        w.writerow(row)

# TreeTime sniffs the delimiter itself and splits quoted fields containing commas
# (e.g. "Ducks, geese and swans"), shifting columns. Give it a minimal 2-column file.
with open(out_path.replace(".csv", "_dates.csv"), "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(["taxa_name", "collection_date"])
    for row in rows:
        w.writerow([row.get("taxa_name", ""), row.get("collection_date", "")])

print(f"Rows with unknown/blank collection_date (will be estimated by TreeTime): {n_unknown}")
PY

# Strip Windows carriage returns (M1/Windows editing artifacts)
perl -i -pe 's/\r//g' "$META_CLEAN" "${META_CLEAN%.csv}_dates.csv"

# 3. Run TreeTime
log "Estimating molecular clock with TreeTime (local-lab + reference subset)..."
treetime --aln "$LOCAL_ALIGNED" \
         --tree "$LOCAL_TREE" \
         --dates "${META_CLEAN%.csv}_dates.csv" \
         --name-column "taxa_name" \
         --date-column "collection_date" \
         --outdir "$TIMETREE_OUT" \
         --confidence \
         --reroot "best"

# 4. Final status check
if [ -d "$TIMETREE_OUT" ]; then
    log "Success: local-lab time-resolved tree generated in $TIMETREE_OUT"
    log "Molecular clock results saved to: $TIMETREE_OUT/molecular_clock.txt"
else
    log "ERROR: TreeTime failed to create output directory."
    exit 1
fi

log "=== Local-lab TreeTime complete for $STAMP ==="
