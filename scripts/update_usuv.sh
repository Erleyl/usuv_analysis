#!/usr/bin/env bash
set -euo pipefail

# ----------------------------
# 0. Setup & Versioning
# ----------------------------
STAMP="${RUN_DATE:-$(date +%Y%m%d)}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
BASE_DIR="$PROJECT_ROOT/usuv_data"
CURATED_DIR="$BASE_DIR/curated"
# Host dictionary (latin -> english name); accept the legacy filename too
BIRD_DICT="$PROJECT_ROOT/dictionary.tsv"
[[ -f "$BIRD_DICT" ]] || BIRD_DICT="$PROJECT_ROOT/bird_dictionary.tsv"
LAST_RUN_FILE="$BASE_DIR/LAST_RUN"

# Filenames for the Master Database
MASTER_FASTA="$CURATED_DIR/usuv_master.fasta"
MASTER_META="$CURATED_DIR/usuv_master_metadata.csv"
LOGFILE="$BASE_DIR/usuv_pipeline_$STAMP.log"

# Load API Key
if [ -f "$PROJECT_ROOT/.env" ]; then export $(grep -v '^#' "$PROJECT_ROOT/.env" | xargs); fi
API_KEY="${NCBI_API_KEY:-}"

mkdir -p "$CURATED_DIR" "$BASE_DIR/raw_runs/$STAMP"
log() { echo "$(date +"[%Y-%m-%d %H:%M:%S]") [FETCH] $*" | tee -a "$LOGFILE"; }

log "=== Starting USUV Incremental Update for $STAMP ==="

# ----------------------------
# 1. Determine Date Range
# ----------------------------
LAST_RUN_DATE="1900-01-01"
[ -f "$LAST_RUN_FILE" ] && LAST_RUN_DATE=$(cat "$LAST_RUN_FILE")
# Force mode (main.sh --force) or manual override: ignore LAST_RUN and re-download
LAST_RUN_DATE="${FETCH_SINCE:-$LAST_RUN_DATE}"

# ----------------------------
# 2. Download New Data
# ----------------------------
RUN_DIR="$BASE_DIR/raw_runs/$STAMP"
log "Downloading genomes released after: $LAST_RUN_DATE"

datasets download virus genome taxon "Usutu virus" \
    --released-after "$LAST_RUN_DATE" \
    --filename "$RUN_DIR/usuv.zip" \
    --api-key "$API_KEY"

# ----------------------------
# 3. Unpack & Path Discovery
# ----------------------------
if [ ! -f "$RUN_DIR/usuv.zip" ]; then log "No usuv.zip found."; exit 0; fi
unzip -o "$RUN_DIR/usuv.zip" -d "$RUN_DIR/ncbi_dataset" >/dev/null
DATADIR=$(find "$RUN_DIR/ncbi_dataset" -type d -name data | head -n 1 || true)
if [ -z "$DATADIR" ]; then log "Data directory not found."; exit 1; fi

# ----------------------------
# 4. Processing with Bird Dictionary & ISO Mapping
# ----------------------------
TEMP_CSV="$RUN_DIR/temp_new.csv"
TEMP_FASTA="$RUN_DIR/temp_new.fasta"

log "Processing raw data with bird dictionary..."

python3 <<PY
import json, pandas as pd
import os, re, sys
from pathlib import Path
from datetime import datetime

# naming rule (country codes, cleaning of date ranges such as "[2015 TO 2019]") lives in sync_master.py
sys.dont_write_bytecode = True      # no __pycache__ folder in scripts/
sys.path.insert(0, "$SCRIPT_DIR")
from sync_master import build_taxa_name

def safe(obj, *keys):
    for k in keys:
        if not isinstance(obj, dict) or k not in obj: return None
        obj = obj[k]
    return obj

bird_map = {}
dict_path = Path("$BIRD_DICT")
if dict_path.exists():
    bird_df = pd.read_csv(dict_path, sep='\t')
    bird_df = bird_df.dropna(subset=['species_latin_name', 'species_name_english'])
    bird_map = {str(k).strip().lower(): str(v).strip()
                for k, v in zip(bird_df['species_latin_name'], bird_df['species_name_english'])}
else:
    print(f"WARNING: host dictionary not found at {dict_path}; species_name_english will be 'unknown'")

meta_dict = {}
jsonl_path = Path("$DATADIR") / "data_report.jsonl"
dl_time = datetime.now().strftime("%Y-%m-%d")

if jsonl_path.exists():
    with open(jsonl_path) as fh:
        for line in fh:
            rec = json.loads(line)
            acc_full = rec.get("accession", "")
            acc_base = acc_full.split('.')[0]
            host_raw = safe(rec, "host", "organismName") or "unknown"
            english_name = bird_map.get(host_raw.strip().lower(), "unknown")
            geo_raw = safe(rec, "location", "geographicLocation") or "unknown"
            coll_date = safe(rec, "isolate", "collectionDate") or "unknown"
            # accession.host.country3.date ; the date is cleaned (a range gives 2015to2019), the
            # collection_date column below keeps the original NCBI value
            taxa_name = build_taxa_name(acc_base, host_raw, geo_raw, coll_date)
            
            # Extract sequence/isolate name safely from nested isolate object
            seq_name = safe(rec, "isolate", "name") or rec.get("description") or rec.get("title") or rec.get("name") or "unknown"
            
            meta_dict[acc_base] = {
                "taxa_name": taxa_name,
                "accession": acc_full,
                "completeness": rec.get("completeness", "unknown"),
                "length": rec.get("length"),
                "organism": safe(rec, "virus", "organismName"),
                "name": seq_name,
                "host": host_raw,
                "species_name_english": english_name,
                "collection_date": coll_date,
                "geo_location": geo_raw,
                "submitter_affiliation": safe(rec, "submitter", "affiliation") or "unknown",
                "submitter_country": safe(rec, "submitter", "country") or "unknown",
                "release_date": rec.get("releaseDate"),
                "update_date": rec.get("updateDate"),
                "download_date": dl_time
            }

with open("$TEMP_FASTA", "w") as fout:
    for fna in Path("$DATADIR").rglob("*.fna"):
        with open(fna) as fin:
            for line in fin:
                if line.startswith(">"):
                    raw_id = line.split()[0][1:]
                    acc_key = re.split(r'[:.]', raw_id)[0]
                    fout.write(f">{meta_dict[acc_key]['taxa_name']}\n") if acc_key in meta_dict else fout.write(f">{acc_key}\n")
                else:
                    fout.write(line)

pd.DataFrame.from_dict(meta_dict, orient='index').to_csv("$TEMP_CSV", index=False)
PY

# ----------------------------
# 5. Incremental Merge & Deduplication
# ----------------------------
log "Merging new data into Master database..."

if [ ! -s "$TEMP_CSV" ]; then
    log "No new data found since last run. Skipping master database merge."
else
    MASTER_DUP_AUDIT="$RUN_DIR/master_duplicate_audit_$STAMP.tsv"

    # Merge FASTA entries
    cat "$TEMP_FASTA" >> "$MASTER_FASTA" 2>/dev/null || cp "$TEMP_FASTA" "$MASTER_FASTA"

    # Python Metadata Merge & Deduplication Step
    python3 <<PY
import pandas as pd
import os
from pathlib import Path

temp_csv = Path("$TEMP_CSV")
master_csv = Path("$MASTER_META")
master_fasta = Path("$MASTER_FASTA")
audit_file = Path("$MASTER_DUP_AUDIT")

# Safely append new metadata using pandas to avoid CSV quoting/delimiter errors
if temp_csv.exists() and temp_csv.stat().st_size > 0:
    new_data = pd.read_csv(temp_csv)
    if not master_csv.exists():
        new_data.to_csv(master_csv, index=False)
    else:
        # Align columns by name (master may predate newer columns such as "name")
        old = pd.read_csv(master_csv, on_bad_lines='skip')
        pd.concat([old, new_data], ignore_index=True, sort=False).to_csv(master_csv, index=False)

if master_csv.exists():
    # on_bad_lines='skip' handles any pre-existing broken lines gracefully
    df = pd.read_csv(master_csv, on_bad_lines='skip')
    
    # Audit duplicate accessions
    duplicates = df[df.duplicated(subset=['accession'], keep=False)]
    if not duplicates.empty:
        duplicates.to_csv(audit_file, sep='\t', index=False)
    
    # Keep the latest record by update_date or download_date
    sort_cols = [col for col in ['update_date', 'download_date'] if col in df.columns]
    if sort_cols:
        df = df.sort_values(by=sort_cols, ascending=False)
    
    df = df.drop_duplicates(subset=['accession'], keep='first')
    df.to_csv(master_csv, index=False)

# Deduplicate FASTA entries cleanly
if master_fasta.exists():
    seen_headers = set()
    unique_records = []
    current_header = None
    current_seq = []
    
    with open(master_fasta, 'r') as fin:
        for line in fin:
            if line.startswith(">"):
                if current_header and current_header not in seen_headers:
                    seen_headers.add(current_header)
                    unique_records.append(f"{current_header}\n{''.join(current_seq)}")
                current_header = line.strip()
                current_seq = []
            else:
                current_seq.append(line)
        if current_header and current_header not in seen_headers:
            seen_headers.add(current_header)
            unique_records.append(f"{current_header}\n{''.join(current_seq)}")
            
    with open(master_fasta, 'w') as fout:
        fout.write("\n".join(unique_records) + "\n" if unique_records else "")
PY
fi

# ----------------------------
# 5b. Repair / sync the master (runs every time, also when there is no new data)
#     - rebuilds taxa_names with unsafe characters (date ranges such as "[2015 TO 2019]")
#     - keeps ONE sequence per accession (the newest) and gives it the metadata taxa_name
#     Old files are kept as *.bak_sync_<time>.
# ----------------------------
if [ -f "$MASTER_FASTA" ] && [ -f "$MASTER_META" ]; then
    log "Syncing master FASTA with master metadata..."
    python3 "$SCRIPT_DIR/sync_master.py" --meta "$MASTER_META" --fasta "$MASTER_FASTA" 2>&1 | tee -a "$LOGFILE"
fi

# ----------------------------
# 6. Final Consistency Verification
# ----------------------------
echo "$(date +"%Y-%m-%d")" > "$LAST_RUN_FILE"
log "=== USUV Incremental Update Completed Successfully ==="