#!/usr/bin/env bash
#
# local_lab_tree.sh
# ==================
# Builds an ML phylogeny using local-lab-curated USUV sequences PLUS a
# curated set of NCBI lineage-reference sequences, so the tree has enough
# context to place local isolates within known USUV lineages.
#
# Rather than re-aligning from scratch, this subsets the alignment already
# produced by align.sh. Sequences are selected in two ways:
#   1. Local lab sequences: matched via data_source=local_lab* in the
#      merged metadata (assigned by merge_lab.sh).
#   2. Reference sequences: matched by ACCESSION against a user-supplied
#      list (reference_lineages.tsv), because taxa_name formatting can
#      drift but the accession is stable. Matching strips NCBI version
#      suffixes (.1, .2, ...) on both sides before comparing.
#
# Prerequisites: run --merge and --align first for the target $RUN_DATE.
#
# Usage:
#   bash local_lab_tree.sh
#   RUN_DATE=20260730 bash local_lab_tree.sh                     # specific run
#   REFERENCE_LIST=/path/to/refs.tsv bash local_lab_tree.sh      # custom ref list
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
ALIGNED_FASTA="$RESULTS_DIR/USUV_Aligned_$STAMP.fasta"

# Reference accession list: accession[TAB]lineage[TAB]year[TAB]location
REFERENCE_LIST="${REFERENCE_LIST:-$PROJECT_ROOT/reference_lineages.tsv}"

COMBINED_IDS="$LOCAL_DIR/local_lab_plus_refs_taxa_ids.txt"
MATCH_REPORT="$LOCAL_DIR/reference_match_report_$STAMP.tsv"
COMBINED_ALIGNED="$LOCAL_DIR/USUV_LocalLab_WithRefs_Aligned_$STAMP.fasta"
PRE="$LOCAL_DIR/USUV_LocalLab_WithRefs_IQTree_$STAMP"

mkdir -p "$LOCAL_DIR"
log() { echo "$(date +"[%Y-%m-%d %H:%M:%S]") [LOCAL-TREE] $*" | tee -a "$BASE_DIR/usuv_pipeline_$STAMP.log"; }

# 1. Pre-flight checks
if [[ ! -f "$MERGED_META" ]]; then
    log "ERROR: Merged metadata not found: $MERGED_META. Run 'main.sh --merge' first."
    exit 1
fi
if [[ ! -f "$ALIGNED_FASTA" ]]; then
    log "ERROR: Aligned FASTA not found: $ALIGNED_FASTA. Run 'main.sh --align' first."
    exit 1
fi
if [[ ! -f "$REFERENCE_LIST" ]]; then
    log "WARNING: Reference list not found at $REFERENCE_LIST -- proceeding with local-lab sequences only, no lineage context."
fi

# 2. Resolve taxa_name for local-lab rows + accession-matched reference rows
log "Selecting local-lab + reference records from $(basename "$MERGED_META")..."
python3 - "$MERGED_META" "$REFERENCE_LIST" "$COMBINED_IDS" "$MATCH_REPORT" <<'PY'
import csv, sys, os

meta_path, ref_path, out_ids, out_report = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]

def base_accession(acc):
    """Strip NCBI version suffix: KC754958.1 -> KC754958"""
    return (acc or "").strip().split(".")[0].upper()

with open(meta_path, newline="", encoding="utf-8") as fh:
    sample = fh.read(4096)
    fh.seek(0)
    delim = ";" if sample.count(";") > sample.count(",") else ","
    rows = list(csv.DictReader(fh, delimiter=delim))

# Local-lab sequences (existing behaviour)
local_ids = [
    row["taxa_name"].strip()
    for row in rows
    if row.get("data_source", "").strip().startswith("local_lab")
    and row.get("taxa_name", "").strip()
]

# Index every metadata row by its version-stripped accession, so a bare
# reference accession (KC754958) matches a versioned one (KC754958.1)
acc_index = {}
for row in rows:
    for key in (row.get("accession", ""), row.get("taxa_name", "").split(".")[0]):
        b = base_accession(key)
        if b:
            acc_index.setdefault(b, row)

ref_ids = []
report_rows = []
if os.path.exists(ref_path):
    with open(ref_path, newline="", encoding="utf-8") as fh:
        ref_rdr = csv.DictReader(fh, delimiter="\t")
        for ref in ref_rdr:
            ref_acc = base_accession(ref.get("accession", ""))
            if not ref_acc:
                continue
            match = acc_index.get(ref_acc)
            if match:
                taxa = match.get("taxa_name", "").strip()
                ref_ids.append(taxa)
                report_rows.append((ref_acc, ref.get("lineage", ""), "FOUND", taxa))
            else:
                report_rows.append((ref_acc, ref.get("lineage", ""), "MISSING", ""))

combined = sorted(set(local_ids) | set(ref_ids))
with open(out_ids, "w", encoding="utf-8") as fh:
    fh.write("\n".join(combined) + ("\n" if combined else ""))

with open(out_report, "w", encoding="utf-8") as fh:
    fh.write("accession\tlineage\tstatus\tmatched_taxa_name\n")
    for r in report_rows:
        fh.write("\t".join(r) + "\n")

n_missing = sum(1 for r in report_rows if r[2] == "MISSING")
print(f"Local-lab sequences: {len(local_ids)}")
print(f"Reference sequences requested: {len(report_rows)}")
print(f"Reference sequences matched: {len(ref_ids)}")
print(f"Reference sequences MISSING (not in master/merged metadata): {n_missing}")
print(f"Combined total: {len(combined)}")
PY

COUNT=$(grep -c . "$COMBINED_IDS" 2>/dev/null || echo 0)
if [[ "$COUNT" -eq 0 ]]; then
    log "ERROR: No local-lab or reference records resolved. Check $MERGED_META and $REFERENCE_LIST."
    exit 1
fi
log "Combined sequences to include: $COUNT"
if [[ -f "$MATCH_REPORT" ]]; then
    MISSING=$(awk -F'\t' 'NR>1 && $3=="MISSING"{print $1}' "$MATCH_REPORT")
    if [[ -n "$MISSING" ]]; then
        log "WARNING: reference accession(s) not found in merged metadata (not yet fetched from NCBI, or version mismatch):"
        while read -r acc; do log "         - $acc"; done <<< "$MISSING"
    fi
fi

# 3. Subset the existing alignment (preserves alignment columns/gaps --
#    no need to re-run MAFFT)
log "Subsetting alignment to local-lab + reference sequences..."
seqkit grep -f "$COMBINED_IDS" "$ALIGNED_FASTA" > "$COMBINED_ALIGNED"

SUBSET_COUNT=$(grep -c ">" "$COMBINED_ALIGNED" || echo 0)
log "Sequences written to subset alignment: $SUBSET_COUNT"
if [[ "$SUBSET_COUNT" -lt "$COUNT" ]]; then
    log "WARNING: $((COUNT - SUBSET_COUNT)) taxa_name(s) resolved above were not found in the aligned FASTA."
    log "         (They may have been dropped by qc_filter.sh, or the header differs from taxa_name.)"
fi
if [[ "$SUBSET_COUNT" -lt 4 ]]; then
    log "ERROR: IQ-TREE needs at least 4 sequences for a meaningful tree; only $SUBSET_COUNT found."
    exit 1
fi

# 4. Build the tree
log "Running IQ-TREE on local-lab + reference alignment (MFP + 1000 UFBoot + 1000 SH-aLRT)..."
REDO=""; [[ "${FORCE:-0}" == "1" ]] && REDO="-redo"
iqtree -s "$COMBINED_ALIGNED" -m MFP -bb 1000 -alrt 1000 -nt AUTO -pre "$PRE" $REDO

log "Done."
log "  Tree file:        ${PRE}.treefile"
log "  Report:            ${PRE}.iqtree"
log "  Consensus:         ${PRE}.contree"
log "  Reference matches: $MATCH_REPORT"
