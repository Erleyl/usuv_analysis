#!/usr/bin/env python3
"""
make_tip_labels.py
==================
Builds a tip-label table from the merged USUV metadata, to load into PearTree (or FigTree,
ggtree, ...) as annotations. Several label styles are made at the same time, so the label
shown on the tips can be changed in the viewer without touching the tree.

The first column is `taxa_name`, the name used for the tips in the pipeline (sequence headers,
alignment, trees), so the table matches the tree. The trees and FASTA files are not changed.

Modelled on the WNV script add_fasta_labels.py (Erley Lizarazo, 2026-08-14), which builds
several label columns from a metadata table.

Columns of the table
  taxa_name, accession, data_source, host, collection_date, geo_location   as in the metadata
  source, project, public, lab_id, country, location, year                 cleaned, for colouring
  label_<style>                                                            one column per style

Fields available in a template
  {acc}          accession without version for NCBI (PX453694.1 -> PX453694); lab IDs unchanged
  {accession}    accession exactly as in the metadata
  {country}      country name (NCBI: first part of geo_location; lab: sheet Country)
  {country3}     3-letter code (nld, deu, ...)
  {location}     detail below the country (NCBI: after the ":" in geo_location; lab: sheet location),
                 "unknown" if there is none
  {host}         host, Latin name as in the metadata
  {host_english} English host name
  {year}         first 4-digit year of the collection date (a range such as [2015 TO 2019] gives 2015)
  {date}         collection date as in the metadata (YYYY-MM-DD, YYYY-MM, YYYY or a cleaned range)
  {source}       whom the sample belongs to: EMC_lab, RIVM, ... or NCBI (not in our sheets).
                 Lab/partner samples that are already public on NCBI count as EMC_lab / RIVM
  {project}      Project of the sheet (DeadWildBirds, TidalParks, rivm, ...), also for the lab's public
                 sequences; NCBI for all other NCBI sequences
  {public}       yes = the sequence is on NCBI, no = only in our lab data
  {lab_id}       the lab's own ID of the sample (empty for other NCBI sequences)
  {name}         isolate name from the metadata
  Missing values become "unknown". Spaces become "_", other characters outside
  A-Z a-z 0-9 . _ - are removed from the values (not from the template text).

Examples
  python3 scripts/make_tip_labels.py --list-styles
  python3 scripts/make_tip_labels.py --metadata usuv_data/curated/20261009/USUV_Merged_20261009_metadata.csv
      -> usuv_data/results/20261009/USUV_tip_labels_20261009.tsv   (same file as main.sh --labels)
  python3 scripts/make_tip_labels.py --template "{source}_{acc}_{country3}_{year}" --lower
"""
import argparse
import csv
import glob
import os
import re
import sys
from collections import Counter

csv.field_size_limit(sys.maxsize)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # project root

# Same codes as merge_lab.sh / update_usuv.sh; other countries: first 3 letters
COUNTRY3 = {
    "netherlands": "nld", "germany": "deu", "belgium": "bel", "france": "fra",
    "austria": "aut", "italy": "ita", "spain": "esp", "united kingdom": "gbr",
    "czech republic": "cze", "hungary": "hun", "poland": "pol", "switzerland": "che",
    "south africa": "zaf", "serbia": "srb", "croatia": "hrv", "slovenia": "svn",
}

STYLES = {
    "acc_country_year":       "{acc}.{country}.{year}",
    "acc_country_host_year":  "{acc}.{country}.{host}.{year}",
    "acc_country_date":       "{acc}_{country}_{date}",
    "acc_host_country3_year": "{acc}.{host}.{country3}.{year}",
    "source_acc_country_year": "{source}.{acc}.{country}.{year}",
}
FIELDS = ["acc", "accession", "country", "country3", "location", "host", "host_english",
          "year", "date", "source", "project", "public", "lab_id", "name"]
UNSAFE = set("()[]:;,'\" \t\r\n")          # would break a label in Newick / FASTA / most viewers
MISSING = {"", "unknown", "nan", "none", "null", "na"}


def clean(value, lower=False):
    v = "" if value is None else str(value).strip()
    if v.lower() in MISSING:
        return "unknown"
    v = re.sub(r"\s+", "_", v)
    v = re.sub(r"[^A-Za-z0-9._-]", "", v)
    if not v:
        return "unknown"
    return v.lower() if lower else v


def load_owners(path):
    """accession (no version) -> owner, for lab/partner samples whose sequence is already on NCBI.
    Read from the selection audit of `main.sh --lab` (columns Source, Project, SeqID_original, GenBank,
    ncbi_exact_accession): the accession listed in the sheet's GenBank column, or an NCBI record with
    the same ID as the lab sample."""
    owners = {}
    with open(path, encoding="utf-8-sig", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            accs = [r.get("GenBank") or ""] + (r.get("ncbi_exact_accession") or "").split(";")
            for a in accs:
                a = a.strip().split(".")[0]
                if a and a not in owners:
                    src, lid = (r.get("Source") or "").strip(), (r.get("SeqID_original") or "").strip()
                    owners[a] = {"source": clean(src) if src else "lab", "project": clean(r.get("Project")),
                                 "lab_id": clean(lid) if lid else ""}
    return owners


def fields_of(row, lower=False, owners=None):
    ds = (row.get("data_source") or "").strip()
    is_ncbi = ds == "ncbi"
    accession = (row.get("accession") or "").strip() or (row.get("taxa_name") or "").split(".")[0]
    acc = re.sub(r"\.\d+$", "", accession) if is_ncbi else accession   # lab IDs such as BRS.24140 keep their dot

    if is_ncbi:
        country_raw = (row.get("geo_location") or "").split(":")[0].strip()
    else:
        country_raw = (row.get("lab_Country") or row.get("submitter_country") or "").strip()
    country = clean(country_raw, lower)
    # location = detail below the country: NCBI "Italy: VR" -> VR; lab rows: the sheet location
    # (merge_lab.sh puts Location/province in geo_location, or the country when there is none)
    geo = (row.get("geo_location") or "").strip()
    if is_ncbi:
        loc_raw = geo.split(":", 1)[1].strip() if ":" in geo else ""
    else:
        loc_raw = "" if geo.lower() == country_raw.lower() else geo
    c3 = COUNTRY3.get(country_raw.lower()) or (re.sub(r"[^a-z]", "", country_raw.lower())[:3]) or "unk"

    date_raw = (row.get("collection_date") or "").strip()
    m = re.search(r"(?:19|20)\d{2}", date_raw)

    # source / project / lab_id say whom the sample belongs to, not where the sequence was downloaded:
    #  - lab/partner rows: Source and Project of the sheet (the merge keeps them as lab_Source, lab_Project)
    #  - NCBI rows that --lab matched to the lab sheet (same ID, or the sheet lists the GenBank accession)
    #    belong to the lab/partner as well, even though they are public
    #  - all other NCBI rows: NCBI
    owner = (owners or {}).get(accession.split(".")[0]) if is_ncbi else None
    if owner:
        source, project, lab_id = owner["source"], owner["project"], owner["lab_id"]
    elif is_ncbi:
        source, project, lab_id = "NCBI", "NCBI", ""
    else:
        src = (row.get("lab_Source") or "").strip()
        source = clean(src) if src and src.lower() != "unknown" else "lab"
        project = clean(row.get("lab_Project"))
        lid = (row.get("lab_SeqID_original") or row.get("lab_SeqID") or "").strip()
        lab_id = clean(lid) if lid and lid.lower() != "unknown" else ""
    public = "yes" if is_ncbi else "no"      # is the sequence on NCBI?

    return {
        "acc": clean(acc), "accession": clean(accession), "country": country, "country3": c3,
        "location": clean(loc_raw, lower),
        "host": clean(row.get("host"), lower), "host_english": clean(row.get("species_name_english"), lower),
        "year": m.group(0) if m else "unknown", "date": clean(date_raw),
        "source": source, "project": project, "public": public, "lab_id": lab_id,
        "name": clean(row.get("name")),
    }


def make_unique(labels):
    seen, out = Counter(), []
    for lab in labels:
        seen[lab] += 1
        out.append(lab if seen[lab] == 1 else f"{lab}_{seen[lab]}")
    return out, sum(1 for n in seen.values() if n > 1)


def main():
    ap = argparse.ArgumentParser(description="Tip-label table for USUV trees (for PearTree and similar viewers).")
    ap.add_argument("--metadata", default=os.path.join(ROOT, "usuv_data", "curated", "USUV_Merged_latest_metadata.csv"),
                    help="merged metadata CSV (needs taxa_name; default: %(default)s)")
    ap.add_argument("--out", help="output TSV. Default: usuv_data/results/<date>/USUV_tip_labels_<date>.tsv "
                                  "when --metadata is USUV_Merged_<date>_metadata.csv (the same file "
                                  "main.sh --labels writes), otherwise usuv_data/results/USUV_tip_labels.tsv")
    ap.add_argument("--template", help='extra style, e.g. "{source}_{acc}_{country3}_{year}" (column label_custom)')
    ap.add_argument("--lab-audit", help="selection audit of main.sh --lab; marks NCBI sequences that belong to the "
                                        "lab/partners. Default: newest usuv_data/lab_sequences/USUV_lab_*_selection_audit.tsv")
    ap.add_argument("--lower", action="store_true", help="lower-case host and country")
    ap.add_argument("--keep-cols", default="", help="extra metadata columns for the table, comma-separated")
    ap.add_argument("--list-styles", action="store_true")
    a = ap.parse_args()

    styles = dict(STYLES)
    if a.template:
        styles["custom"] = a.template
    if a.list_styles:
        print("Styles:")
        for k, v in styles.items():
            print(f"  {k:26s} {v}")
        print("Fields:", ", ".join("{" + f + "}" for f in FIELDS))
        return 0

    if not os.path.exists(a.metadata):
        sys.exit(f"Metadata not found: {a.metadata}")
    with open(a.metadata, encoding="utf-8-sig", newline="") as fh:
        rows = list(csv.DictReader(fh))
    if not rows or "taxa_name" not in rows[0]:
        sys.exit("The metadata needs a taxa_name column.")

    # whom do the public (NCBI) sequences belong to?
    audit = a.lab_audit
    if not audit:
        found = sorted(glob.glob(os.path.join(ROOT, "usuv_data", "lab_sequences", "USUV_lab_*_selection_audit.tsv")))
        audit = found[-1] if found else None
    owners = load_owners(audit) if audit and os.path.exists(audit) else {}
    if owners:
        print(f"lab audit {os.path.basename(audit)}: {len(owners)} public accessions belong to lab/partner samples")
    else:
        print("WARNING: no lab selection audit found (run main.sh --lab): public sequences of the lab show as NCBI")

    # build every style
    labels = {k: [] for k in styles}
    plain = {"source": [], "project": [], "public": [], "lab_id": [], "country": [], "location": [], "year": []}
    for r in rows:
        f = fields_of(r, a.lower, owners)
        for k in plain:
            plain[k].append(f[k])
        for k, tpl in styles.items():
            try:
                labels[k].append(tpl.format(**f))
            except KeyError as e:
                sys.exit(f"Style '{k}': unknown field {e}. Available: {', '.join(FIELDS)}")
    n_dup = {}
    for k in styles:
        labels[k], n_dup[k] = make_unique(labels[k])
        bad = [x for x in labels[k] if UNSAFE & set(x)]
        if bad:
            sys.exit(f"Style '{k}' gives labels with characters that break viewers/Newick "
                     f"(any of ( ) [ ] : ; , ' \" or spaces), e.g. {bad[0]!r}")

    out = a.out
    if not out:
        m = re.match(r"^USUV_Merged_(.+)_metadata\.csv$", os.path.basename(a.metadata))
        if m and m.group(1) != "latest":
            out = os.path.join(ROOT, "usuv_data", "results", m.group(1), f"USUV_tip_labels_{m.group(1)}.tsv")
        else:
            out = os.path.join(ROOT, "usuv_data", "results", "USUV_tip_labels.tsv")
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    keep = ["taxa_name", "accession", "data_source", "host", "collection_date", "geo_location"] + \
           [c for c in a.keep_cols.split(",") if c.strip()]
    keep = [c for c in dict.fromkeys(keep) if c in rows[0]]
    with open(out, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(keep + list(plain) + [f"label_{k}" for k in styles])
        for i, r in enumerate(rows):
            w.writerow([r.get(c, "") for c in keep] + [plain[k][i] for k in plain] + [labels[k][i] for k in styles])
    print(f"{len(rows)} rows, {len(styles)} label styles -> {out}")
    for k in styles:
        if n_dup[k]:
            print(f"  note: style 'label_{k}' had {n_dup[k]} repeated labels; later ones got _2, _3, ...")
    return 0


if __name__ == "__main__":
    sys.exit(main())
