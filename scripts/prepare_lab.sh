#!/usr/bin/env bash
# ==============================================================================
# Script Name: prepare_lab.sh            (main.sh step: --lab)
# Description: Prepares the local-lab and partner sequences for the merge step.
#              Runs scripts/prepare_lab_fasta.R (clean sequences, give every row a
#              safe unique ID, drop sequences that are already on NCBI, add partner
#              sheets such as RIVM) and installs the result where merge_lab.sh reads
#              its "raw" lab inputs:
#                  usuv_data/local_lab/internal_genomes.fasta
#                  usuv_data/local_lab/internal_metadata.csv
#              merge_lab.sh then builds the NCBI-style names (SeqID.host.country.year)
#              and renames the FASTA headers itself, as in earlier runs.
#              NOT used on purpose: usuv_data/lab_*_curated.* (from curate_lab_data.sh).
#              merge_lab.sh prefers those files and in that mode expects FASTA headers
#              that are already NCBI-style names, which this step does not produce.
#              The full set of outputs (incl. the "_all" export and the selection
#              audit) is kept in usuv_data/lab_sequences/USUV_lab_<RUN_DATE>*.
#
#              Length / N content are NOT filtered here: that is done once, for lab
#              and NCBI sequences together, by qc_filter.sh (main.sh --qc).
#
# Settings (environment variables, all optional):
#   LAB_SHEET     lab overview CSV. Default: newest USUVSequencingOverview_ED*.csv
#                 in usuv_data/local_lab
#   LAB_EXTERNAL  partner sheets as FILE=LABEL[,FILE=LABEL]. Default:
#                 usuv_data/local_lab/rivm_seqs.csv=RIVM (if that file exists).
#                 Use LAB_EXTERNAL=none for no partner data.
#   LAB_R_ARGS    extra options for prepare_lab_fasta.R, for example
#                 LAB_R_ARGS="--emma-filter routine --emma-extra 2025"
#
# Credits: the study selection and metadata logic come from Emmanuelle Munger's
#          R scripts (scripts/emma/); see the header of prepare_lab_fasta.R.
# ==============================================================================
set -euo pipefail

STAMP="${RUN_DATE:-$(date +%Y%m%d)}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

BASE_DIR="$PROJECT_ROOT/usuv_data"
LAB_DIR="$BASE_DIR/local_lab"
OUT_DIR="$BASE_DIR/lab_sequences"
CURATED_DIR="$BASE_DIR/curated"

# NCBI reference (from the --fetch step), used to find lab sequences that are already public
NCBI_META="$CURATED_DIR/usuv_master_metadata.csv"
NCBI_FASTA="$CURATED_DIR/usuv_master.fasta"

# What merge_lab.sh reads in "raw" mode (canonical names, preferred over any other file in local_lab)
RAW_LAB_FASTA="$LAB_DIR/internal_genomes.fasta"
RAW_LAB_META="$LAB_DIR/internal_metadata.csv"
# merge_lab.sh would prefer these over the files above, so they must not exist
CURATED_LAB_FASTA="$BASE_DIR/lab_sequences_curated.fasta"
CURATED_LAB_META="$BASE_DIR/lab_metadata_curated.csv"

mkdir -p "$OUT_DIR" "$BASE_DIR/logs"
log()  { echo "$(date +"[%Y-%m-%d %H:%M:%S]") [LAB] $*" | tee -a "$BASE_DIR/usuv_pipeline_$STAMP.log"; }
fail() { log "ERROR: $*"; exit 1; }

command -v Rscript >/dev/null 2>&1 || fail "Rscript not found (install R with the packages dplyr and stringr)."
Rscript -e 'suppressPackageStartupMessages({library(dplyr); library(stringr)})' >/dev/null 2>&1 \
    || fail "R packages missing. In R: install.packages(c('dplyr','stringr'))"

# 1. Inputs ---------------------------------------------------------------------
if [[ -n "${LAB_SHEET:-}" ]]; then
    SHEET="$LAB_SHEET"
else
    # ED<yymmdd> sorts chronologically, so the last one is the newest sheet
    SHEET="$(ls -1 "$LAB_DIR"/USUVSequencingOverview_ED*.csv 2>/dev/null | sort | tail -n 1 || true)"
fi
[[ -n "$SHEET" && -f "$SHEET" ]] || fail "No lab overview CSV found (looked for $LAB_DIR/USUVSequencingOverview_ED*.csv). Set LAB_SHEET=/path/to/sheet.csv"
log "Lab sheet: $SHEET"

if [[ -z "${LAB_EXTERNAL+x}" ]]; then
    [[ -f "$LAB_DIR/rivm_seqs.csv" ]] && LAB_EXTERNAL="$LAB_DIR/rivm_seqs.csv=RIVM" || LAB_EXTERNAL="none"
fi

R_ARGS=(--input "$SHEET" --outdir "$OUT_DIR" --date-tag "$STAMP")

if [[ "$LAB_EXTERNAL" != "none" && -n "$LAB_EXTERNAL" ]]; then
    R_ARGS+=(--external "$LAB_EXTERNAL")
    log "Partner sheet(s): $LAB_EXTERNAL"
else
    log "No partner sheets."
fi

if [[ -f "$NCBI_META" ]]; then
    R_ARGS+=(--ncbi-metadata "$NCBI_META")
else
    log "WARNING: $NCBI_META not found (run --fetch first). Lab sequences that are already on NCBI will NOT be detected by ID."
fi
if [[ -f "$NCBI_FASTA" ]]; then
    R_ARGS+=(--ncbi-fasta "$NCBI_FASTA")
else
    log "WARNING: $NCBI_FASTA not found. Identical-sequence check against NCBI skipped."
fi

# 2. Run the R script -------------------------------------------------------------
log "Running prepare_lab_fasta.R ..."
# LAB_R_ARGS is split on spaces on purpose
# shellcheck disable=SC2086
Rscript "$SCRIPT_DIR/prepare_lab_fasta.R" "${R_ARGS[@]}" ${LAB_R_ARGS:-} 2>&1 | tee -a "$BASE_DIR/usuv_pipeline_$STAMP.log"

PREFIX="$OUT_DIR/USUV_lab_$STAMP"
[[ -s "$PREFIX.fasta" && -s "${PREFIX}_metadata.csv" ]] || fail "prepare_lab_fasta.R did not produce $PREFIX.fasta / _metadata.csv"

# 3. Install for merge_lab.sh (keep the previous version as .prev) ---------------
if [[ -f "$CURATED_LAB_FASTA" || -f "$CURATED_LAB_META" ]]; then
    fail "Found $CURATED_LAB_FASTA / $CURATED_LAB_META (from curate_lab_data.sh). merge_lab.sh would use them instead of the files prepared here. Move them away first, e.g.: mkdir -p $BASE_DIR/old_curated && mv $BASE_DIR/lab_*_curated.* $BASE_DIR/old_curated/"
fi
for pair in "$PREFIX.fasta|$RAW_LAB_FASTA" "${PREFIX}_metadata.csv|$RAW_LAB_META"; do
    src="${pair%%|*}"; dst="${pair##*|}"
    if [[ -f "$dst" ]] && ! cmp -s "$src" "$dst"; then
        cp -p "$dst" "$dst.prev"
        log "Previous $(basename "$dst") saved as $(basename "$dst").prev"
    fi
    cp "$src" "$dst"
done

N_SEQ="$(grep -c '^>' "$RAW_LAB_FASTA" || true)"
log "Installed $N_SEQ lab/partner sequences for the merge:"
log "  $RAW_LAB_FASTA"
log "  $RAW_LAB_META"
log "Full outputs and the per-row audit: $OUT_DIR/USUV_lab_${STAMP}*"
log "Next: bash main.sh --merge --qc  (merge_lab.sh uses the files above automatically)"
