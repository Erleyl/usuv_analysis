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
BIRD_DICT="$PROJECT_ROOT/bird_dictionary.tsv"
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
import os, re
from pathlib import Path
from datetime import datetime

ISO_MAP = {"Italy": "ita", "Germany": "deu", "Austria": "aut", "Hungary": "hun", "Spain": "esp", "France": "fra", "Netherlands": "nld", "Belgium": "bel", "Switzerland": "che", "Czech Republic": "cze", "Poland": "pol", "USA": "usa", "South Africa": "zaf", "Serbia": "srb", "Croatia": "hrv"}

def safe(obj, *keys):
    for k in keys:
        if not isinstance(obj, dict) or k not in obj: return None
        obj = obj[k]
    return obj

bird_map = {}
dict_path = Path("$BIRD_DICT")
if dict_path.exists():
    bird_df = pd.read_csv(dict_path, sep='\t')
    bird_map = dict(zip(bird_df['species_latin_name'], bird_df['species_name_english']))

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
            english_name = bird_map.get(host_raw, "unknown")
            geo_raw = safe(rec, "location", "geographicLocation") or "unknown"
            coll_date = safe(rec, "isolate", "collectionDate") or "unknown"
            country_name = geo_raw.split(":")[0].strip() if ":" in geo_raw else geo_raw
            country3 = ISO_MAP.get(country_name, country_name[:3].lower())
            host_clean = re.sub(r'[^a-z0-9_]', '', host_raw.lower().replace(" ", "_"))
            taxa_name = f"{acc_base}.{host_clean}.{country3}.{coll_date}"
            
            meta_dict[acc_base] = {
                "taxa_name": taxa_name, "accession": acc_full, "completeness": rec.get("completeness", "unknown"),
                "length": rec.get("length"), "organism": safe(rec, "virus", "organismName"),
                "host": host_raw, "species_name_english": english_name, 
                "collection_date": coll_date, "geo_location": geo_raw,
                "submitter_affiliation": safe(rec, "submitter", "affiliation") or "unknown",
                "submitter_country": safe(rec, "submitter", "country") or "unknown",
                "release_date": rec.get("releaseDate"), "update_date": rec.get("updateDate"), "download_date": dl_time
            }

with open("$TEMP_FASTA", "w") as fout:
    for fna in Path("$DATADIR").rglob("*.fna"):
        with open(fna) as fin:
            for line in fin:
                if line.startswith(">"):
                    raw_id = line.split()[0][1:]
                    acc_key = re.split(r'[:.]', raw_id)[0]
                    fout.write(f">{meta_dict[acc_key]['taxa_name']}\n") if acc_key in meta_dict else fout.write(f">{acc_key}\n")
                else: fout.write(line)

pd.DataFrame.from_dict(meta_dict, orient='index').to_csv("$TEMP_CSV", index=False)
PY

# ----------------------------
# 5. Incremental Merge & Deduplication
# ----------------------------
log "Merging new data into Master database..."

# Merge Metadata
if [ ! -f "$MASTER_META" ]; then
    cp "$TEMP_CSV" "$MASTER_META"
else
    tail -n +2 "$TEMP_CSV" >> "$MASTER_META"
fi

# Merge FASTA
cat "$TEMP_FASTA" >> "$MASTER_FASTA" 2>/dev/null || cp "$TEMP_FASTA" "$MASTER_FASTA"

# Python Deduplication Step
python3 <<PY
import pandas as pd
import os
# Read and remove duplicates based on taxa_name, keeping the most recent download
df = pd.read_csv("$MASTER_META").drop_duplicates(subset=["taxa_name"], keep="last")
df.to_csv("$MASTER_META", index=False)

valid_taxa = set(df["taxa_name"].tolist())
seen_taxa = set()

with open("${MASTER_FASTA}.tmp", "w") as fout, open("$MASTER_FASTA", "r") as fin:
    valid = False
    for line in fin:
        if line.startswith(">"):
            h = line.strip()[1:]
            valid = h in valid_taxa and h not in seen_taxa
            if valid:
                seen_taxa.add(h)
                fout.write(line)
        elif valid:
            fout.write(line)
os.replace("${MASTER_FASTA}.tmp", "$MASTER_FASTA")
PY

# ----------------------------
# 6. Final Consistency Verification
# ----------------------------
echo "$(date +"%Y-%m-%d")" > "$LAST_RUN_FILE"

log "--- Final Verification ---"
META_COUNT=$(tail -n +2 "$MASTER_META" 2>/dev/null | wc -l | xargs || echo "0")
FASTA_COUNT=$(grep -c "^>" "$MASTER_FASTA" 2>/dev/null | xargs || echo "0")

log "Consistency Check: Master Metadata records ($META_COUNT) vs Master FASTA sequences ($FASTA_COUNT)"

if [ "$META_COUNT" -eq "$FASTA_COUNT" ]; then
    log "Success: Master Database is synchronized."
else
    log "Warning: Metadata/FASTA count mismatch. Check logs."
fi

log "=== Update Complete for $STAMP ==="