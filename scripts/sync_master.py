#!/usr/bin/env python3
"""
sync_master.py
==============
Keeps the NCBI master database consistent (used by update_usuv.sh, step --fetch).

Problem it solves
  NCBI sometimes gives a collection date as a RANGE, e.g. "[2015 TO 2019]". The sequence name
  taxa_name = accession.host.country.date then contained a space and brackets: the space cut the
  FASTA header, the bracket starts a comment in Newick, and the merge could no longer match the
  sequence to its metadata row (it created an "unknown" placeholder instead). In addition, when an
  NCBI record changes between downloads (date, host name, ...) its name changes, and the master
  FASTA, which was de-duplicated by header, kept the old and the new copy of the same accession.

What it does (safe to run any number of times)
  1. metadata: a taxa_name that contains characters other than A-Z a-z 0-9 . _ - is rebuilt from the
     row (host, geo_location, collection_date). Names that are fine are never changed.
     The collection_date column keeps the original value.
  2. FASTA: one record per accession (the last one = the newest download), and its header is set to
     the taxa_name of the metadata row, so FASTA and metadata always match.
  3. Before anything is rewritten, the old files are copied to <file>.bak_sync_<time>.

Naming rule (also imported by update_usuv.sh)
  accession(without version).host(lower, a-z0-9_).country3.date
  date: YYYY-MM-DD / YYYY-MM / YYYY unchanged; a range "[2015 TO 2019]" -> 2015to2019;
        anything else is reduced to A-Za-z0-9_- ; nothing left -> unknown.

Usage
  python3 scripts/sync_master.py                      # repairs usuv_data/curated/usuv_master*
  python3 scripts/sync_master.py --dry-run            # only report
  python3 scripts/sync_master.py --meta M.csv --fasta F.fasta
"""
import argparse
import csv
import os
import re
import shutil
import sys
from datetime import datetime

csv.field_size_limit(sys.maxsize)

ISO_MAP = {"Italy": "ita", "Germany": "deu", "Austria": "aut", "Hungary": "hun", "Spain": "esp",
           "France": "fra", "Netherlands": "nld", "Belgium": "bel", "Switzerland": "che",
           "Czech Republic": "cze", "Poland": "pol", "USA": "usa", "South Africa": "zaf",
           "Serbia": "srb", "Croatia": "hrv"}
UNSAFE = re.compile(r"[^A-Za-z0-9._-]")
RANGE = re.compile(r"^\[?\s*(\d{4}(?:-\d{2}(?:-\d{2})?)?)\s+TO\s+(\d{4}(?:-\d{2}(?:-\d{2})?)?)\s*\]?$", re.I)


def clean_date_for_name(value):
    v = "" if value is None else str(value).strip()
    if v == "" or v.lower() in ("unknown", "nan", "none", "null"):
        return "unknown"
    m = RANGE.match(v)
    if m:
        return f"{m.group(1)}to{m.group(2)}"
    v = re.sub(r"\s+", "_", v)
    return re.sub(r"[^A-Za-z0-9_-]", "", v) or "unknown"


def build_taxa_name(acc_base, host_raw, geo_raw, coll_date):
    host_clean = re.sub(r"[^a-z0-9_]", "", str(host_raw).lower().replace(" ", "_"))
    geo_raw = str(geo_raw)
    country_name = geo_raw.split(":")[0].strip() if ":" in geo_raw else geo_raw
    country3 = re.sub(r"[^a-z]", "", ISO_MAP.get(country_name, country_name[:3].lower())) or "unk"
    return f"{acc_base}.{host_clean}.{country3}.{clean_date_for_name(coll_date)}"


def read_fasta(path):
    """list of (header_without_>, sequence_lines_text)"""
    recs, hdr, seq = [], None, []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if hdr is not None:
                    recs.append((hdr, "".join(seq)))
                hdr, seq = line[1:].strip(), []
            elif hdr is not None:
                seq.append(line if line.endswith("\n") else line + "\n")
    if hdr is not None:
        recs.append((hdr, "".join(seq)))
    return recs


def sync(meta_path, fasta_path, dry_run=False):
    with open(meta_path, encoding="utf-8-sig", newline="") as fh:
        rdr = csv.DictReader(fh)
        fields = rdr.fieldnames
        rows = list(rdr)

    # 1. metadata: rebuild only the names that contain unsafe characters
    renamed_meta = 0
    for r in rows:
        name = (r.get("taxa_name") or "").strip()
        if UNSAFE.search(name) or not name:
            acc = (r.get("accession") or "").split(".")[0]
            r["taxa_name"] = build_taxa_name(acc, r.get("host") or "unknown",
                                             r.get("geo_location") or "unknown", r.get("collection_date"))
            renamed_meta += 1
    name_of = {(r.get("accession") or "").split(".")[0]: r["taxa_name"] for r in rows}

    # 2. FASTA: one record per accession (the last one), header = metadata taxa_name
    recs = read_fasta(fasta_path)
    keep, renamed_fa, dropped = {}, 0, 0
    for i, (hdr, seq) in enumerate(recs):
        key = re.split(r"[:.]", hdr.split()[0])[0]
        if key not in name_of:
            key = hdr.split()[0]                       # not in the metadata: dedupe by header only
        if key in keep:
            dropped += 1
        keep[key] = (i, hdr, seq)
    out = []
    for key, (i, hdr, seq) in sorted(keep.items(), key=lambda kv: kv[1][0]):
        new = name_of.get(key)
        if new and hdr.split()[0] != new:
            renamed_fa += 1
            hdr = new
        out.append((hdr, seq))
    missing = sum(1 for k in name_of if k not in keep)

    print(f"[sync] metadata: {len(rows)} rows, {renamed_meta} taxa_name rebuilt (unsafe characters)")
    print(f"[sync] FASTA: {len(recs)} records -> {len(out)}; {dropped} older copies of the same accession removed, "
          f"{renamed_fa} headers set to the metadata name; {missing} metadata rows without a sequence")
    if renamed_meta == 0 and dropped == 0 and renamed_fa == 0:
        print("[sync] nothing to change")
        return 0
    if dry_run:
        print("[sync] dry run: no file written")
        return 0

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    for p in (meta_path, fasta_path):
        shutil.copy2(p, f"{p}.bak_sync_{stamp}")
    tmp_m, tmp_f = meta_path + ".tmp", fasta_path + ".tmp"
    with open(tmp_m, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    with open(tmp_f, "w", encoding="utf-8") as fh:
        for hdr, seq in out:
            fh.write(f">{hdr}\n{seq}")
    os.replace(tmp_m, meta_path)
    os.replace(tmp_f, fasta_path)
    print(f"[sync] written; previous files kept as *.bak_sync_{stamp}")
    return 0


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    cur = os.path.join(root, "usuv_data", "curated")
    ap = argparse.ArgumentParser(description="Repair the USUV master database (names, duplicate sequences).")
    ap.add_argument("--meta", default=os.path.join(cur, "usuv_master_metadata.csv"))
    ap.add_argument("--fasta", default=os.path.join(cur, "usuv_master.fasta"))
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    for p in (a.meta, a.fasta):
        if not os.path.exists(p):
            print(f"[sync] {p} not found, nothing to do")
            return 0
    return sync(a.meta, a.fasta, a.dry_run)


if __name__ == "__main__":
    sys.exit(main())
