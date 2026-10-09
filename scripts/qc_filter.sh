#!/usr/bin/env bash
set -euo pipefail

STAMP="${RUN_DATE:-$(date +%Y%m%d)}"
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/usuv_data"
INPUT="$BASE_DIR/curated/$STAMP/USUV_Merged_$STAMP.fasta"
# Legacy flat layout fallback (runs created before dated folders)
[[ -f "$INPUT" ]] || INPUT="$BASE_DIR/curated/USUV_Merged_$STAMP.fasta"
[[ -f "$INPUT" ]] || { echo "[QC] ERROR: merged FASTA for $STAMP not found. Run --merge first." >&2; exit 1; }
OUTPUT="$BASE_DIR/results/$STAMP/USUV_QC_Passed_$STAMP.fasta"
mkdir -p "$(dirname "$OUTPUT")"

# Re-apply exclusion lists so nothing excluded reaches alignment/tree steps,
# even if the merged FASTA was built before a list was updated.
PROJECT_ROOT="$(dirname "$BASE_DIR")"
FILTERED_INPUT="$(dirname "$OUTPUT")/USUV_Merged_excl_$STAMP.fasta"
INPUT_FASTA="$INPUT" OUT_FASTA="$FILTERED_INPUT" \
EXC_ACC="$PROJECT_ROOT/exclude_accessions.txt" EXC_SEQS="$PROJECT_ROOT/exclude_seqs.txt" python3 - <<'PY'
import os

def load(path, loose):
    out = []
    if os.path.exists(path):
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                line = line.split("#")[0].strip() if not loose else line.strip()
                if line:
                    out.append(line if loose else line.split()[0])
    return out

acc = load(os.environ["EXC_ACC"], True)    # same loose rules as merge_lab.sh
seqs = set(load(os.environ["EXC_SEQS"], False))  # exact IDs

def drop(h):
    if h in seqs or h.split(".")[0] in seqs:
        return True
    h_clean = h.split("|")[0].strip()
    return any(h_clean == t or t in h or h_clean in t for t in acc)

kept = removed = 0
write = False
with open(os.environ["INPUT_FASTA"], encoding="utf-8") as fin, open(os.environ["OUT_FASTA"], "w", encoding="utf-8") as fout:
    for line in fin:
        if line.startswith(">"):
            write = not drop(line[1:].strip().split()[0])
            kept += write
            removed += (not write)
        if write:
            fout.write(line)
print(f"[QC] Exclusion lists re-applied: removed {removed}, kept {kept}")
PY
INPUT="$FILTERED_INPUT"

# Default thresholds follow Emmanuelle Munger's rule (240729_PrepSeqData.R):
# keep sequences longer than 8800 nt with fewer than 2200 N, i.e. about 20% N
# for an 11 kb genome. (Before: MIN_LEN=10000, MAX_NS=10.0.)
# This is the ONLY N/length filter for lab and NCBI sequences: the lab curation
# step (main.sh --lab: prepare_lab_fasta.R) does not filter on length or N.
# Override per run: MIN_LEN=10000 MAX_NS=10 ./main.sh --qc
MIN_LEN="${MIN_LEN:-8800}"
MAX_NS="${MAX_NS:-20.0}"
QC_AUDIT="$(dirname "$OUTPUT")/USUV_QC_dropped_$STAMP.tsv"

echo "[QC] Filtering: Min Length ${MIN_LEN}bp, Max N content ${MAX_NS}% (letter N only)"
IN_FASTA="$INPUT" OUT_FASTA="$OUTPUT" AUDIT="$QC_AUDIT" MIN_LEN="$MIN_LEN" MAX_NS="$MAX_NS" python3 - <<'PY'
import os

min_len = int(float(os.environ["MIN_LEN"]))
max_ns = float(os.environ["MAX_NS"])

def records(path):
    hdr, chunks = None, []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if hdr is not None:
                    yield hdr, chunks
                hdr, chunks = line.rstrip("\n"), []
            elif hdr is not None:
                chunks.append(line)
    if hdr is not None:
        yield hdr, chunks

kept = dropped_len = dropped_n = 0
with open(os.environ["OUT_FASTA"], "w", encoding="utf-8") as fout, \
     open(os.environ["AUDIT"], "w", encoding="utf-8") as audit:
    audit.write("id\tlength\tn_count\tn_percent\treason\n")
    for hdr, chunks in records(os.environ["IN_FASTA"]):
        seq = "".join(c.strip() for c in chunks).upper()
        length = len(seq)
        n_count = seq.count("N")
        n_pct = 100.0 * n_count / length if length else 100.0
        sid = hdr[1:].split()[0]
        reasons = []
        if length < min_len:
            reasons.append(f"length<{min_len}")
        if n_pct > max_ns:
            reasons.append(f"N>{max_ns:g}%")
        if reasons:
            audit.write(f"{sid}\t{length}\t{n_count}\t{n_pct:.2f}\t{';'.join(reasons)}\n")
            dropped_len += reasons[0].startswith("length")
            dropped_n += any(r.startswith("N>") for r in reasons)
            continue
        kept += 1
        fout.write(hdr + "\n" + "".join(chunks))

print(f"[QC] Kept {kept}; dropped for length: {dropped_len}; dropped for N content: {dropped_n} (overlap possible)")
PY

echo "[QC] Done. Sequences remaining: $(grep -c ">" "$OUTPUT")"
echo "[QC] Dropped-sequence audit: $QC_AUDIT"
