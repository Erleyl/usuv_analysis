#!/usr/bin/env bash
# ==============================================================================
# Script Name: merge_lab.sh
# Date: 2026-06-11
# Author: Erley Lizarazo
# Description: Merges master NCBI sequence databases with local lab datasets.
#              Systematically drops identical genomic replicates via precise
#              header token and substring mapping rules.
# ==============================================================================
set -euo pipefail

# 1. Environment & Path Setup
STAMP="${RUN_DATE:-$(date +%Y%m%d)}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

BASE_DIR="$PROJECT_ROOT/usuv_data"
CURATED_DIR="$BASE_DIR/curated"
EXCLUDE_FILE="$PROJECT_ROOT/exclude_accessions.txt"
# Local-lab sequences already on NCBI (exact-ID match, lab sequences only)
EXCLUDE_SEQS_FILE="$PROJECT_ROOT/exclude_seqs.txt"
[[ -f "$EXCLUDE_SEQS_FILE" ]] || : > "$EXCLUDE_SEQS_FILE"

# Master files from fetch step and new dated snapshots
MASTER_FASTA="$CURATED_DIR/usuv_master.fasta"
MASTER_META="$CURATED_DIR/usuv_master_metadata.csv"
# Dated snapshot folder: curated/<STAMP>/ (master files stay in curated/)
RUN_CURATED_DIR="$CURATED_DIR/$STAMP"
FINAL_FASTA="$RUN_CURATED_DIR/USUV_Merged_$STAMP.fasta"
FINAL_META="$RUN_CURATED_DIR/USUV_Merged_${STAMP}_metadata.csv"
LATEST_META="$CURATED_DIR/USUV_Merged_latest_metadata.csv"
LAB_HEADER_MAP="$RUN_CURATED_DIR/USUV_Merged_${STAMP}_lab_header_map.tsv"
DUP_AUDIT="$RUN_CURATED_DIR/USUV_Merged_${STAMP}_duplicate_sequences_audit.tsv"

# Preferred curated lab files produced by curate_lab_data.sh integration
# (old route, now in scripts/old_scripts/; main.sh --lab installs the raw files below instead)
CURATED_LAB_FASTA="$BASE_DIR/lab_sequences_curated.fasta"
CURATED_LAB_META="$BASE_DIR/lab_metadata_curated.csv"

# Fallback raw lab files
RAW_LAB_DIR="$BASE_DIR/local_lab"
RAW_LAB_FASTA="$RAW_LAB_DIR/internal_genomes.fasta"
RAW_LAB_META="$RAW_LAB_DIR/internal_metadata.csv"

LAB_FASTA=""
LAB_META=""
LAB_MODE="none"
MISSING_METADATA_POLICY="${MISSING_METADATA_POLICY:-placeholder}" # placeholder|drop

mkdir -p "$RAW_LAB_DIR" "$RUN_CURATED_DIR"
log() { echo "$(date +"[%Y-%m-%d %H:%M:%S]") [MERGE] $*" | tee -a "$BASE_DIR/usuv_pipeline_$STAMP.log"; }

resolve_lab_inputs() {
    # Prefer curated, NCBI-compatible inputs if present
    if [[ -f "$CURATED_LAB_FASTA" && -f "$CURATED_LAB_META" ]]; then
        LAB_FASTA="$CURATED_LAB_FASTA"
        LAB_META="$CURATED_LAB_META"
        LAB_MODE="curated"
        log "Using curated lab inputs (NCBI-compatible)."
        return
    fi

    # Fallback to raw local_lab files
    if [[ -f "$RAW_LAB_FASTA" ]]; then
        LAB_FASTA="$RAW_LAB_FASTA"
    else
        local first_fasta
        first_fasta=$(find "$RAW_LAB_DIR" -maxdepth 1 -type f \( -name "*.fasta" -o -name "*.fa" -o -name "*.fna" \) | head -n 1 || true)
        if [[ -n "$first_fasta" ]]; then
            LAB_FASTA="$first_fasta"
        fi
    fi

    if [[ -f "$RAW_LAB_META" ]]; then
        LAB_META="$RAW_LAB_META"
    else
        # Try first CSV in local_lab if the canonical name is not present
        local first_csv
        first_csv=$(find "$RAW_LAB_DIR" -maxdepth 1 -type f -name "*.csv" | head -n 1 || true)
        if [[ -n "$first_csv" ]]; then
            LAB_META="$first_csv"
        fi
    fi

    if [[ -n "$LAB_FASTA" && -n "$LAB_META" ]]; then
        LAB_MODE="raw"
        log "Using raw lab inputs from local_lab (will normalize metadata)."
    elif [[ -n "$LAB_FASTA" ]]; then
        LAB_MODE="fasta_only"
        log "Found local FASTA without metadata. Only FASTA will be merged."
    else
        LAB_MODE="none"
        log "No local lab inputs detected. Using NCBI-only dataset."
    fi
}

merge_metadata() {
    log "Merging metadata..."

    if [[ ! -f "$MASTER_META" ]]; then
        log "ERROR: Master metadata not found at $MASTER_META. Run --fetch first."
        exit 1
    fi

    LAB_META_PATH="$LAB_META" LAB_MODE="$LAB_MODE" LAB_FASTA_PATH="$LAB_FASTA" python3 <<PY
import os
import re
import csv
from datetime import datetime

master_meta = r"$MASTER_META"
out_meta = r"$FINAL_META"
latest_meta = r"$LATEST_META"
header_map = r"$LAB_HEADER_MAP"
exclude_file = r"$EXCLUDE_FILE"
exclude_seqs_file = r"$EXCLUDE_SEQS_FILE"
lab_meta_path = os.environ.get("LAB_META_PATH", "")
lab_mode = os.environ.get("LAB_MODE", "none")
lab_fasta_path = os.environ.get("LAB_FASTA_PATH", "")

core_cols = [
    "taxa_name", "accession", "completeness", "length", "organism", "name", "host",
    "species_name_english", "collection_date", "geo_location",
    "submitter_affiliation", "submitter_country", "release_date", "update_date", "download_date"
]

def read_csv(path):
    with open(path, "r", encoding="utf-8", newline="") as fh:
        sample = fh.read(4096)
        fh.seek(0)
        delimiter = ","
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",;\t")
            delimiter = dialect.delimiter
        except Exception:
            delimiter = ";" if sample.count(";") > sample.count(",") else ","
        rdr = csv.DictReader(fh, delimiter=delimiter)
        return list(rdr), (rdr.fieldnames or [])

def write_csv(path, rows, fieldnames):
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)

master_rows, master_fields = read_csv(master_meta)
for row in master_rows:
    for c in core_cols:
        if c not in row or row[c] in (None, ""):
            row[c] = "unknown"
    row["data_source"] = "ncbi"

excluded = []
if os.path.exists(exclude_file):
    with open(exclude_file, "r", encoding="utf-8") as fh:
        excluded = [line.strip() for line in fh if line.strip()]

def is_missing(v):
    if v is None:
        return True
    s = str(v).strip()
    return s == "" or s.lower() == "nan"

def pick(row, keys):
    for k in keys:
        if k in row and not is_missing(row.get(k)):
            return str(row.get(k)).strip()
    return ""

def norm_host(value):
    if is_missing(value):
        return "unknown"
    return re.sub(r"[^a-z0-9_]+", "", str(value).strip().lower().replace(" ", "_")) or "unknown"

def year_from_date(value):
    if is_missing(value):
        return "unknown"
    txt = str(value).strip()
    m = re.search(r"(19|20)\\d{2}", txt)
    return m.group(0) if m else "unknown"

country_map = {
    "netherlands": "nld", "germany": "deu", "belgium": "bel", "france": "fra",
    "austria": "aut", "italy": "ita", "spain": "esp", "united kingdom": "gbr",
    "czech republic": "cze", "hungary": "hun", "poland": "pol", "switzerland": "che",
    "south africa": "zaf", "serbia": "srb", "croatia": "hrv", "slovenia": "svn"
}

def country3(value):
    if is_missing(value):
        return "unk"
    txt = str(value).strip().lower()
    if txt in country_map:
        return country_map[txt]
    return re.sub(r"[^a-z]", "", txt)[:3] or "unk"

def is_excluded_value(value):
    if is_missing(value):
        return False
    s = str(value).strip()
    # Normalize structural delimiters for accurate substring containment checks
    s_clean = s.split()[0].split("|")[0].strip()
    for token in excluded:
        t_clean = token.strip()
        if s_clean == t_clean or t_clean in s or s_clean in t_clean:
            return True
    return False

def fasta_id_index(path):
    if not path or not os.path.exists(path):
        return set()
    out = set()
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                token = line[1:].strip().split()[0]
                out.add(token)
                out.add(token.split(".")[0])
                out.add(re.split(r"[|:.]", token)[0])
    return out

def load_ids(path):
    ids = set()
    if os.path.exists(path):
        with open(path, "r", encoding="utf-8") as fh:
            for line in fh:
                line = line.split("#")[0].strip()
                if line:
                    ids.add(line.split()[0])
    return ids

exclude_seqs = load_ids(exclude_seqs_file)
dropped_lab_dups = 0

def in_exclude_seqs(*names):
    for n in names:
        if n and (n in exclude_seqs or n.split(".")[0] in exclude_seqs):
            return True
    return False

lab_rows = []
lab_fields = []
header_map_rows = []

if lab_mode in {"curated", "raw"} and lab_meta_path and os.path.exists(lab_meta_path):
    raw_rows, raw_fields = read_csv(lab_meta_path)
    lab_fields = raw_fields
    local_fasta_ids = fasta_id_index(lab_fasta_path)

    for idx, row in enumerate(raw_rows, start=1):
        accession = pick(row, ["accession", "SeqID", "seqid", "SequenceID", "id", "ID"])
        taxa_name = pick(row, ["taxa_name"])

        host = pick(row, ["host", "Species_Latin", "Species_English"]) or "unknown"
        species_english = pick(row, ["species_name_english", "Species_English"]) or "unknown"
        collection_date = pick(row, ["collection_date", "DeathOrSampleDate", "AnalysisDate"]) or "unknown"
        country = pick(row, ["submitter_country", "Country"]) or "unknown"
        location = pick(row, ["geo_location", "Location", "province"]) or ""

        if not taxa_name:
            seq_token = accession if accession else f"LAB{idx:06d}"
            taxa_name = f"{seq_token}.{norm_host(host)}.{country3(country)}.{year_from_date(collection_date)}"

        if not accession:
            accession = taxa_name.split(".")[0]

        if lab_mode == "raw" and local_fasta_ids:
            acc_base = accession.split(".")[0]
            if accession not in local_fasta_ids and acc_base not in local_fasta_ids and taxa_name not in local_fasta_ids:
                continue

        out = {c: "unknown" for c in core_cols}
        out["taxa_name"] = taxa_name
        out["accession"] = accession
        out["completeness"] = pick(row, ["completeness"]) or "COMPLETE"
        length_txt = pick(row, ["length", "SEQ_LENGTH"])
        if length_txt:
            try:
                out["length"] = str(int(float(length_txt)))
            except Exception:
                out["length"] = "unknown"
        out["organism"] = pick(row, ["organism"]) or "Usutu virus"
        out["name"] = pick(row, ["name", "isolate", "isolate_name", "seq_name"]) or "unknown"
        out["host"] = host
        out["species_name_english"] = species_english
        out["collection_date"] = collection_date
        out["geo_location"] = location if location else country
        out["submitter_affiliation"] = pick(row, ["submitter_affiliation", "Applicant"]) or "Local Laboratory"
        out["submitter_country"] = country
        out["release_date"] = pick(row, ["release_date"]) or datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")
        out["update_date"] = pick(row, ["update_date"]) or datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")
        out["download_date"] = pick(row, ["download_date"]) or datetime.utcnow().strftime("%Y-%m-%d")
        out["data_source"] = "local_lab"

        for col in raw_fields:
            if col in out:
                continue
            out[f"lab_{col}"] = row.get(col, "")

        if in_exclude_seqs(out["accession"], out["taxa_name"], pick(row, ["SeqID", "seqid"])):
            dropped_lab_dups += 1
            continue

        if not is_excluded_value(out["accession"]) and not is_excluded_value(out["taxa_name"]):
            lab_rows.append(out)
            header_map_rows.append((out["accession"], out["taxa_name"]))

merged = master_rows + lab_rows

# Apply exclusions to merged rows
if excluded:
    kept = []
    for row in merged:
        acc = row.get("accession", "")
        taxa = row.get("taxa_name", "")
        if not (is_excluded_value(acc) or is_excluded_value(taxa)):
            kept.append(row)
    merged = kept

# Deduplicate by taxa_name (keep last)
seen = set()
dedup_reversed = []
for row in reversed(merged):
    key = row.get("taxa_name", "")
    if key not in seen:
        seen.add(key)
        dedup_reversed.append(row)
merged = list(reversed(dedup_reversed))

all_fields = []
for c in core_cols + ["data_source"] + master_fields:
    if c and c not in all_fields:
        all_fields.append(c)
for row in merged:
    for k in row.keys():
        if k and k not in all_fields:
            all_fields.append(k)

write_csv(out_meta, merged, all_fields)
write_csv(latest_meta, merged, all_fields)

with open(header_map, "w", encoding="utf-8", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["source_id", "taxa_name"])
    for src, taxa in header_map_rows:
        w.writerow([src, taxa])

print(f"Merged metadata saved: {out_meta}")
print(f"Latest metadata snapshot saved: {latest_meta}")
print(f"Merged metadata records: {len(merged)}")
print(f"Local-lab rows dropped via exclude_seqs.txt: {dropped_lab_dups}")
PY
}

merge_sequences() {
    log "Merging FASTA sequences..."

    # Start with master FASTA
    cat "$MASTER_FASTA" > "$FINAL_FASTA"

    # Append local FASTA if available
    if [[ -n "$LAB_FASTA" && -f "$LAB_FASTA" ]]; then
        if [[ "$LAB_MODE" == "raw" && -f "$LAB_HEADER_MAP" ]]; then
            local normalized_lab_fasta
            normalized_lab_fasta="$RUN_CURATED_DIR/USUV_Merged_${STAMP}_local_lab_normalized.fasta"

            LAB_FASTA_PATH="$LAB_FASTA" LAB_MAP_PATH="$LAB_HEADER_MAP" LAB_OUT_FASTA="$normalized_lab_fasta" python3 <<PY
import os
import re
import csv

in_fa = os.environ["LAB_FASTA_PATH"]
map_tsv = os.environ["LAB_MAP_PATH"]
out_fa = os.environ["LAB_OUT_FASTA"]

mapping = {}
with open(map_tsv, "r", encoding="utf-8") as fh:
    rdr = csv.DictReader(fh, delimiter="\t")
    for row in rdr:
        src = (row.get("source_id") or "").strip()
        taxa = (row.get("taxa_name") or "").strip()
        if src and taxa:
            mapping[src] = taxa
            mapping[src.split(".")[0]] = taxa

def normalized_id(raw_header: str) -> str:
    token = raw_header.split()[0].strip()
    token_base = re.split(r"[|:.]", token)[0]
    if token in mapping:
        return mapping[token]
    if token_base in mapping:
        return mapping[token_base]
    return token

with open(in_fa, "r", encoding="utf-8") as fin, open(out_fa, "w", encoding="utf-8") as fout:
    for line in fin:
        if line.startswith(">"):
            raw = line[1:].strip()
            fout.write(">" + normalized_id(raw) + "\n")
        else:
            fout.write(line)
PY

            cat "$normalized_lab_fasta" >> "$FINAL_FASTA"
            log "Appended normalized local lab FASTA: $normalized_lab_fasta"
        else
            cat "$LAB_FASTA" >> "$FINAL_FASTA"
            log "Appended local lab FASTA: $LAB_FASTA"
        fi
    fi

    # Apply exclusion list
    if [ -f "$EXCLUDE_FILE" ]; then
        log "Checking exclusion list..."
        COUNT=$(grep -c '^' "$EXCLUDE_FILE" || echo 0)
        if [ "$COUNT" -gt 0 ]; then
            log "Removing $COUNT blacklisted accessions from this run..."
            EXC_FILE="$EXCLUDE_FILE" IN_FASTA="$FINAL_FASTA" OUT_FASTA="$FINAL_FASTA.tmp" python3 <<PY
import os

exc_file = os.environ["EXC_FILE"]
in_fa = os.environ["IN_FASTA"]
out_fa = os.environ["OUT_FASTA"]

with open(exc_file, "r", encoding="utf-8") as fh:
    excluded = [line.strip() for line in fh if line.strip()]

def keep(header: str) -> bool:
    # Isolate sample ID prefix completely before any pipe annotations or whitespace modifications
    h_clean = header.strip().split()[0].split("|")[0].strip()
    for token in excluded:
        t_clean = token.strip()
        if h_clean == t_clean or t_clean in header or h_clean in t_clean:
            return False
    return True

write_seq = False
with open(in_fa, "r", encoding="utf-8") as fin, open(out_fa, "w", encoding="utf-8") as fout:
    for line in fin:
        if line.startswith(">"):
            hdr = line[1:].strip().split()[0]
            write_seq = keep(hdr)
            if write_seq:
                fout.write(line)
        elif write_seq:
            fout.write(line)
PY
            mv "$FINAL_FASTA.tmp" "$FINAL_FASTA"
        fi
    else
        log "No exclusion list found at $EXCLUDE_FILE. Skipping blacklist step."
    fi

    # Drop local-lab sequences already on NCBI (exact ID match; see exclude_seqs.txt)
    if [ -s "$EXCLUDE_SEQS_FILE" ]; then
        EXC_FILE="$EXCLUDE_SEQS_FILE" IN_FASTA="$FINAL_FASTA" OUT_FASTA="$FINAL_FASTA.tmp" python3 <<PY
import os

exc = set()
with open(os.environ["EXC_FILE"], "r", encoding="utf-8") as fh:
    for line in fh:
        line = line.split("#")[0].strip()
        if line:
            exc.add(line.split()[0])

removed = 0
write_seq = False
with open(os.environ["IN_FASTA"], "r", encoding="utf-8") as fin, open(os.environ["OUT_FASTA"], "w", encoding="utf-8") as fout:
    for line in fin:
        if line.startswith(">"):
            hdr = line[1:].strip().split()[0]
            write_seq = not (hdr in exc or hdr.split(".")[0] in exc)
            if not write_seq:
                removed += 1
        if write_seq:
            fout.write(line)
os.replace(os.environ["OUT_FASTA"], os.environ["IN_FASTA"])
print(f"Sequences removed via exclude_seqs.txt: {removed}")
PY
    fi

    log "Final sequence snapshot created: $FINAL_FASTA"
}

deduplicate_fasta_headers() {
    log "De-duplicating FASTA headers (IQ-TREE requires unique names)..."

    IN_FASTA="$FINAL_FASTA" OUT_FASTA="$FINAL_FASTA.tmp" AUDIT_FILE="$DUP_AUDIT" python3 <<PY
import os
import hashlib
from collections import defaultdict

in_fa = os.environ["IN_FASTA"]
out_fa = os.environ["OUT_FASTA"]
audit_file = os.environ["AUDIT_FILE"]

seen = set()
kept = 0
dropped = 0
write_seq = False

# Build duplicate audit first
entries = []
current_header = None
current_seq_chunks = []

with open(in_fa, "r", encoding="utf-8") as fin:
    for line in fin:
        if line.startswith(">"):
            if current_header is not None:
                seq = "".join(current_seq_chunks)
                entries.append((current_header, seq))
            current_header = line[1:].strip().split()[0]
            current_seq_chunks = []
        else:
            current_seq_chunks.append(line.strip())
    if current_header is not None:
        seq = "".join(current_seq_chunks)
        entries.append((current_header, seq))

positions = defaultdict(list)
hashes_by_header = defaultdict(list)
seq_hash_positions = defaultdict(list)
seq_hash_headers = defaultdict(list)
for idx, (hdr, seq) in enumerate(entries, start=1):
    seq_hash = hashlib.sha256(seq.encode("utf-8")).hexdigest()
    positions[hdr].append(idx)
    hashes_by_header[hdr].append(seq_hash)
    seq_hash_positions[seq_hash].append(idx)
    seq_hash_headers[seq_hash].append(hdr)

with open(audit_file, "w", encoding="utf-8") as af:
    af.write(
        "audit_type\tgroup_key\toccurrences\tunique_headers\tunique_sequence_variants\tidentical_sequences\t"
        "all_positions\tall_headers\tkept_header\tkept_position\tdropped_headers\tdropped_positions\n"
    )

    duplicate_header_rows = 0
    for hdr, pos in sorted(positions.items()):
        if len(pos) <= 1:
            continue
        hs = hashes_by_header[hdr]
        uniq = len(set(hs))
        identical = "yes" if uniq == 1 else "no"
        duplicate_header_rows += 1
        kept_pos = pos[0]
        dropped_headers = [hdr for _ in pos[1:]]
        dropped_pos = ",".join(str(p) for p in pos[1:])
        af.write(
            f"duplicate_header\t{hdr}\t{len(pos)}\t1\t{uniq}\t{identical}\t"
            f"{','.join(str(p) for p in pos)}\t{','.join([hdr for _ in pos])}\t{hdr}\t{kept_pos}\t"
            f"{','.join(dropped_headers)}\t{dropped_pos}\n"
        )

    duplicate_sequence_rows = 0
    for seq_hash, pos in sorted(seq_hash_positions.items()):
        if len(pos) <= 1:
            continue
        headers = seq_hash_headers[seq_hash]
        unique_headers = sorted(set(headers))
        kept_pos = pos[0]
        kept_header = headers[0]
        dropped_pos = ",".join(str(p) for p in pos[1:])
        dropped_headers = ",".join(headers[1:])
        duplicate_sequence_rows += 1
        af.write(
            f"duplicate_sequence\tsha256:{seq_hash}\t{len(pos)}\t{len(unique_headers)}\t1\tyes\t"
            f"{','.join(str(p) for p in pos)}\t{','.join(headers)}\t{kept_header}\t{kept_pos}\t"
            f"{dropped_headers}\t{dropped_pos}\n"
        )

with open(in_fa, "r", encoding="utf-8") as fin, open(out_fa, "w", encoding="utf-8") as fout:
    for line in fin:
        if line.startswith(">"):
            hdr = line[1:].strip().split()[0]
            if hdr in seen:
                write_seq = False
                dropped += 1
            else:
                seen.add(hdr)
                write_seq = True
                kept += 1
                fout.write(line)
        elif write_seq:
            fout.write(line)

os.replace(out_fa, in_fa)
print(f"Unique FASTA headers kept: {kept}")
print(f"Duplicate FASTA headers removed: {dropped}")
print(f"Duplicate sequence audit written: {audit_file}")
print(f"Duplicate header groups audited: {duplicate_header_rows}")
print(f"Duplicate sequence-content groups audited: {duplicate_sequence_rows}")
PY
}

verify_merge() {
    local fasta_count
    fasta_count=$(grep -c "^>" "$FINAL_FASTA" 2>/dev/null || echo 0)
    local meta_count
    meta_count=$(tail -n +2 "$FINAL_META" 2>/dev/null | wc -l | xargs || echo 0)

    log "Verification: merged FASTA sequences=$fasta_count, merged metadata rows=$meta_count"
}

sync_metadata_to_fasta() {
    log "Synchronizing merged metadata with merged FASTA headers (policy: $MISSING_METADATA_POLICY)..."

    FINAL_FASTA_PATH="$FINAL_FASTA" FINAL_META_PATH="$FINAL_META" LATEST_META_PATH="$LATEST_META" POLICY="$MISSING_METADATA_POLICY" python3 <<PY
import os
import csv

fasta_path = os.environ["FINAL_FASTA_PATH"]
meta_path = os.environ["FINAL_META_PATH"]
latest_path = os.environ["LATEST_META_PATH"]
policy = os.environ.get("POLICY", "placeholder").strip().lower()

fasta_ids = []
with open(fasta_path, "r", encoding="utf-8") as fh:
    for line in fh:
        if line.startswith(">"):
            fasta_ids.append(line[1:].strip().split()[0])

fasta_set = set(fasta_ids)

with open(meta_path, "r", encoding="utf-8", newline="") as fh:
    sample = fh.read(4096)
    fh.seek(0)
    delimiter = ";" if sample.count(";") > sample.count(",") else ","
    rdr = csv.DictReader(fh, delimiter=delimiter)
    rows = list(rdr)
    fields = rdr.fieldnames or []

filtered = [r for r in rows if (r.get("taxa_name", "").strip() in fasta_set)]

meta_taxa = set(r.get("taxa_name", "").strip() for r in filtered if r.get("taxa_name", "").strip())
missing_from_meta = [h for h in fasta_ids if h not in meta_taxa]

if policy == "placeholder" and missing_from_meta:
    defaults = {
        "completeness": "unknown",
        "length": "unknown",
        "organism": "Usutu virus",
        "host": "unknown",
        "species_name_english": "unknown",
        "collection_date": "unknown",
        "geo_location": "unknown",
        "submitter_affiliation": "unknown",
        "submitter_country": "unknown",
        "release_date": "unknown",
        "update_date": "unknown",
        "download_date": "unknown",
        "data_source": "local_lab_placeholder",
    }
    for header in missing_from_meta:
        accession = header.split(".")[0]
        row = {k: "unknown" for k in fields}
        row["taxa_name"] = header
        if "accession" in row:
            row["accession"] = accession
        for k, v in defaults.items():
            if k in row:
                row[k] = v
        filtered.append(row)
    print(f"Policy placeholder: added {len(missing_from_meta)} metadata placeholder rows")

if policy == "drop" and missing_from_meta:
    keep = set(meta_taxa)
    out_fasta = fasta_path + ".tmp"
    write_seq = False
    kept = 0
    with open(fasta_path, "r", encoding="utf-8") as fin, open(out_fasta, "w", encoding="utf-8") as fout:
        for line in fin:
            if line.startswith(">"):
                hdr = line[1:].strip().split()[0]
                write_seq = hdr in keep
                if write_seq:
                    kept += 1
                    fout.write(line)
            elif write_seq:
                fout.write(line)
    os.replace(out_fasta, fasta_path)
    print(f"Policy drop: removed {len(missing_from_meta)} FASTA headers without metadata (kept {kept})")

for out in (meta_path, latest_path):
    with open(out, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in filtered:
            w.writerow(r)

print(f"Metadata rows before sync: {len(rows)}")
print(f"Metadata rows after sync: {len(filtered)}")
print(f"Missing FASTA headers before policy handling: {len(missing_from_meta)}")
PY
}

# 2. Start building the snapshot
log "Combining Master database with local lab data for version $STAMP..."

if [ ! -f "$MASTER_FASTA" ]; then
    log "ERROR: Master database not found at $MASTER_FASTA. Run --fetch first."
    exit 1
fi

resolve_lab_inputs
merge_metadata
merge_sequences
deduplicate_fasta_headers
sync_metadata_to_fasta
verify_merge

log "Merge completed."
log "- Sequence output: $FINAL_FASTA"
log "- Metadata output: $FINAL_META"
log "- Duplicate audit: $DUP_AUDIT"