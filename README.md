# USUV Virus Phylogenetic Analysis Pipeline

## Description
Data analysis of USUV sequences. Intended to generate ML phylogenetic trees for recurrent data overviews and reports.

The pipeline combines the public NCBI sequences with the sequences of the local lab and of
partners (currently RIVM), removes what is already public, filters on quality, and builds
an alignment, a maximum-likelihood tree and a time tree.

## Pipeline at a glance

```
--fetch   NCBI download (incremental)              -> usuv_data/curated/usuv_master*.{fasta,csv}
--lab     lab sheet + partner sheets (R)           -> usuv_data/local_lab/internal_{genomes.fasta,metadata.csv}
--merge   NCBI master + lab/partner sequences      -> usuv_data/curated/<date>/USUV_Merged_<date>.{fasta,_metadata.csv}
--qc      length / N filter                        -> usuv_data/results/<date>/USUV_QC_Passed_<date>.fasta
--align   MAFFT alignment
--tree    IQ-TREE maximum-likelihood tree
--time    TreeTime molecular clock
--labels  tip-label table (TSV) for PearTree (made automatically by --merge)
```

Optional: `--local-tree` / `--local-time` build a tree of the lab samples plus the reference
lineages listed in `reference_lineages.tsv`.

## Usage

```bash
cd usuv_phylo_analysis
bash main.sh                  # prints the operating manual

# full run, in this order
bash main.sh --fetch --lab --merge --qc --align --tree --time

# lab or partner data changed, no new NCBI download needed
bash main.sh --lab --merge --qc --align --tree

# everything (--lab runs right after --fetch)
bash main.sh --all
```

| Command | Description |
| :--- | :--- |
| `bash main.sh --dry-run` | Checks that the tools (incl. `Rscript`) are installed. |
| `bash main.sh --fetch` | Downloads and curates NCBI data (incremental updates). |
| `bash main.sh --lab` | Prepares lab and partner sequences for the merge (see below). |
| `bash main.sh --merge` | Merges the NCBI master with the lab data and applies the exclusion lists. Also writes the tip-label table (see below). |
| `bash main.sh --qc` | Filters on length and N content (default: at least 8800 bp, at most 20 % N). |
| `bash main.sh --align` | Aligns sequences using MAFFT. |
| `bash main.sh --tree` | Generates the phylogeny with IQ-TREE. |
| `bash main.sh --time` | Runs TreeTime. |
| `bash main.sh --labels` | Remakes the tip-label table (TSV) for PearTree, for example with another `LABEL_TEMPLATE` (see below). |
| `bash main.sh --all` | Runs the whole pipeline (`--lab` after `--fetch`). |
| `--force` | Redo steps that were already run for this `RUN_DATE`. Combine with other options. |

Steps run in the order you give them. Use `RUN_DATE=YYYYMMDD` to work on an older run.

## Lab and partner data (`--lab`)

### Inputs (in `usuv_data/local_lab/`)
- **Lab sheet:** the newest `USUVSequencingOverview_ED<yymmdd>.csv`. It has one row per sequence
  with metadata columns and a `Sequence` column.
- **Partner sheet:** `rivm_seqs.csv` (same columns as the lab sheet, at least `SeqID` and
  `Sequence`). It is added as source `RIVM`; its IDs get the prefix `RIVM_` (`4940` -> `RIVM_4940`).
- **NCBI reference:** `usuv_data/curated/usuv_master_metadata.csv` and `usuv_master.fasta`
  (from `--fetch`), used to find lab sequences that are already public.

### What it does (`scripts/prepare_lab.sh` -> `scripts/prepare_lab_fasta.R`)
1. Cleans the sequences (removes gaps and whitespace, upper case). Text such as `FAILED` in the
   `Sequence` column is not sequence data.
2. Gives every row a safe, unique ID: unsafe characters are removed, a missing `SeqID` is replaced
   by the `UniqueID` (or `LABROW<n>`), repeated IDs get `_dup2`, `_repeat`, ...
3. Detects sequences that are already on NCBI: the lab `SeqID` equals the last `/` part of the NCBI
   `name` (for example `USUV/NL/2017/MEA155` -> `MEA155`).
4. Writes two exports and an audit (see Outputs) and installs the selection for the merge step.

A row goes to the merge unless it
- has a GenBank accession in the sheet, or is already on NCBI;
- is tagged `Excluded_*`, `Exlcuded_*` or `Not80Perc` in `Included_In`;
- belongs to the projects `Culture` or `MosquitoSurveillance_CMV_GenExternal`;
- is the worse of a sample and its repeat (`Remark` = `Repeat`; the sequence with less N is kept,
  and if the first attempt failed the repeat is kept);
- is an exact duplicate (same ID, same sequence) of another row.

Emmanuelle Munger's study selection (`Included_In` list) is **not** applied by default. Switch it on with
`LAB_R_ARGS="--emma-filter routine --emma-extra 2025"`.

Length and N content are **not** filtered here, `--qc` does that once for all sequences.

### Outputs
In `usuv_data/lab_sequences/`, with prefix `USUV_lab_<date>`:

| File | Content |
| :--- | :--- |
| `_all.fasta`, `_all_metadata.csv` | Every row with real sequence data, nothing filtered out, with flag columns. |
| `.fasta`, `_metadata.csv` | The selection sent to the merge (same IDs as above). |
| `_ENA_submission_metadata.csv` | ENA-style table for own lab rows (partner data is not included). |
| `_selection_audit.tsv` | One line per input row: what happened and why. |

The selection is installed as `usuv_data/local_lab/internal_genomes.fasta` and `internal_metadata.csv`
(the previous version is kept as `.prev`). `merge_lab.sh` reads these files and builds the sequence
names (`SeqID.host.country.year`, for example
`AS201600289_BEC009.hippoboscidae_strix_nebulosa.nld.2016`) from the sheet columns.

Review columns in the metadata: `seqid_note`, `coordinates_note`, `n_percent`,
`ncbi_identical_seq` (sequence identical to an NCBI record), `identical_sequence_as`
(identical to another row), `ncbi_base_id_match_REVIEW` (same ID after removing a `_suffix`).
Confidential coordinates are rounded to 2 decimals; impossible coordinates (for example
`5.19848e+16`, a lost decimal separator) are removed and the sequence is kept.

### Settings (environment variables)
| Variable | Default |
| :--- | :--- |
| `LAB_SHEET` | newest `USUVSequencingOverview_ED*.csv` in `usuv_data/local_lab/` |
| `LAB_EXTERNAL` | `usuv_data/local_lab/rivm_seqs.csv=RIVM` if that file exists; `none` = no partner data; several: `FILE=LABEL,FILE=LABEL` |
| `LAB_R_ARGS` | none. Extra options for `prepare_lab_fasta.R` (see the header of that script) |

Example: `LAB_EXTERNAL=none LAB_R_ARGS="--keep-repeats" bash main.sh --lab`

## Quality filter (`--qc`)
`scripts/qc_filter.sh` re-applies the exclusion lists and then keeps sequences with at least
`MIN_LEN` bp and at most `MAX_NS` percent N. The defaults follow Emmanuelle Munger's rule
(more than 8800 nt and fewer than 2200 N, about 20 % for an 11 kb genome):

```bash
MIN_LEN=8800 MAX_NS=20 bash main.sh --qc      # the defaults
MIN_LEN=10000 MAX_NS=10 bash main.sh --qc     # the previous, stricter values
```

Everything that is dropped is listed with the reason in `usuv_data/results/<date>/USUV_QC_dropped_<date>.tsv`.
The thresholds apply to NCBI and lab sequences alike.

## Tip labels for PearTree (made by `--merge`)
Inside the pipeline the sequence headers, metadata and trees are joined on `taxa_name`
(`accession.host.country.date`), so the trees keep that name. For figures, a table with readable
labels is made from the merged metadata and loaded into PearTree next to the tree. The trees and FASTA
files are not changed. The table is made automatically at the end of every `--merge` (it only needs the
merged metadata, not the tree); if that fails, the pipeline just continues with a warning. `--labels` remakes it.

```bash
bash main.sh --labels
LABEL_TEMPLATE="{source}_{acc}_{country3}_{year}" LABEL_LOWER=1 bash main.sh --labels   # extra own style
```

It writes `usuv_data/results/<date>/USUV_tip_labels_<date>.tsv`, one row per sequence of the merged data:

| Column(s) | Content |
| :--- | :--- |
| `taxa_name` | first column, equal to the tip names of the trees |
| `accession`, `data_source`, `host`, `collection_date`, `geo_location` | as in the merged metadata |
| `source`, `project`, `public`, `lab_id` | whom the sample belongs to: `source` = `EMC_lab`, `RIVM` or `NCBI` (not in our sheets), `project` = Project of the sheet (`DeadWildBirds`, `TidalParks`, `rivm`, ...; `NCBI` for other public sequences), `lab_id` = our own ID. Lab/partner samples that are already public on NCBI keep their source, project and `lab_id`; `public` = `yes` if the sequence is on NCBI, `no` if it is only in our lab data |
| `country`, `location`, `year` | cleaned values, handy for colouring tips (`location` = detail below the country, `unknown` if none) |
| `label_acc_country_year`, `label_acc_country_host_year`, `label_acc_country_date`, `label_acc_host_country3_year`, `label_source_acc_country_year` (+ `label_custom`) | ready-made tip labels in several styles |

Which sequences belong to the lab or to a partner is taken from the selection audit of `--lab` (`usuv_data/lab_sequences/USUV_lab_<date>_selection_audit.tsv`): the accession in the sheet's `GenBank` column, or an NCBI record with the same ID as the lab sample. Run `--lab` before `--labels`; without that audit the lab's public sequences show as `NCBI`.

Load the TSV as annotations in PearTree (the first column is matched with the tip names), then choose a
`label_*` column for the tip text and, for example, `project`, `source`, `public`, `country`, `location`, `year` or `host`
for colours. Labels are unique and contain only letters, digits and `. _ -`.

The script can also be used directly (`scripts/make_tip_labels.py`, same modelling as the WNV script
`add_fasta_labels.py`):

```bash
python3 scripts/make_tip_labels.py --list-styles          # styles and template fields
python3 scripts/make_tip_labels.py --metadata usuv_data/curated/<date>/USUV_Merged_<date>_metadata.csv   # -> usuv_data/results/<date>/USUV_tip_labels_<date>.tsv (same file as --labels)
```

## Master database repair (`sync_master.py`)
`--fetch` ends with `scripts/sync_master.py`, which keeps `usuv_master.fasta` and `usuv_master_metadata.csv`
consistent. It runs every time (also when there is nothing new to download) and can be run alone:

```bash
python3 scripts/sync_master.py --dry-run     # only report
python3 scripts/sync_master.py               # repair usuv_data/curated/usuv_master*
```

- **Date ranges.** NCBI sometimes gives a collection date as a range (`[2015 TO 2019]`). The sequence name
  `accession.host.country.date` then used to contain a space and brackets, the merge could not match the
  sequence to its metadata and created an `unknown` placeholder. Now the range is written as `2015to2019` in
  the name; the `collection_date` column keeps the original value. Other names are never changed.
- **One sequence per accession.** When an NCBI record changes between downloads (date, host name, ...) its
  name changes, and the master used to keep the old and the new copy. Now only the newest copy is kept and its
  header is always the `taxa_name` of the metadata row.
- **TreeTime.** `--time` and `--local-time` write date ranges as blank (unknown) in their dates file, because
  TreeTime cannot read the `[A TO B]` syntax.
- The old files are kept as `usuv_master*.bak_sync_<time>`. Merged results made before this fix
  (`usuv_data/curated/<date>/`) still contain the old names until `--fetch` and then `--merge --force` are run again.

## Exclusion lists
- `exclude_accessions.txt`: NCBI accessions (or lab IDs) that must never reach the alignment.
- `exclude_seqs.txt`: exact IDs of lab sequences that are already on NCBI.

Both are applied in `--merge` and again in `--qc`.

## Installation

```bash
conda env create -f environment.yaml
conda activate usuv_env
```

The environment includes R (`r-base`, `r-dplyr`, `r-stringr`), which is needed for `--lab`.
TreeTime is not in the environment file and has to be installed separately.

## Directory structure

```
usuv_phylo_analysis/
├── main.sh                      # entry point
├── environment.yaml
├── dictionary.tsv               # host names (Dutch / Latin / English)
├── exclude_accessions.txt
├── exclude_seqs.txt
├── reference_lineages.tsv       # references for --local-tree
├── scripts/
│   ├── update_usuv.sh           # --fetch
│   ├── sync_master.py           # --fetch: repairs names / duplicate sequences in the master
│   ├── prepare_lab.sh           # --lab (wrapper)
│   ├── prepare_lab_fasta.R      # --lab (cleaning, selection, NCBI duplicate check)
│   ├── merge_lab.sh             # --merge
│   ├── qc_filter.sh             # --qc
│   ├── align.sh                 # --align
│   ├── phylogeny.sh             # --tree
│   ├── timetree.sh              # --time
│   ├── local_lab_tree.sh        # --local-tree
│   ├── local_lab_timetree.sh    # --local-time
│   ├── tip_labels.sh            # --labels (optional): table for PearTree
│   ├── make_tip_labels.py       # used by tip_labels.sh, also usable directly
│   ├── emma/                    # original R scripts by Emmanuelle Munger (reference only)
│   └── old_scripts/             # retired scripts, see the README in that folder
└── usuv_data/
    ├── local_lab/               # lab sheets, partner sheets, installed internal_* files
    ├── lab_sequences/           # outputs of --lab (exports, audit)
    ├── curated/                 # NCBI master + dated merged FASTA / metadata
    ├── results/<date>/          # QC result, alignment, trees
    └── logs/
```

## Known issues (accepted for now)

- **Date detail differs between sources:** NCBI names carry the full date, lab names only the year.
- **Host names of lab sequences** come from `Species_Latin` in the sheet. Where that cell holds an English name
  or `Unknown`, the host part of the name is not the Latin name.
- **`KC754958`** (Africa 1) is in `exclude_accessions.txt` and in `reference_lineages.tsv`, so `--local-tree`
  always reports it as MISSING.
- **Impossible coordinates in the lab sheet** (for example `5.19848e+16`, a lost decimal separator) are removed by
  `--lab` (column `coordinates_note`) while the sequence is kept. They have to be corrected in the sheet.

## Troubleshooting

1. **`Rscript not found` / `R packages missing` at `--lab`**: install R and the packages
   (`conda env update -f environment.yaml`, or `install.packages(c("dplyr","stringr"))` in R).
2. **`--lab` stops with "Found .../lab_sequences_curated.fasta"**: these files come from the old
   curation route. `merge_lab.sh` would use them instead of the files prepared by `--lab`, so move
   them away (the error message shows the command).
3. **A lab sequence is missing from the tree**: look it up in
   `usuv_data/lab_sequences/USUV_lab_<date>_selection_audit.tsv` (reason for each row) and in
   `usuv_data/results/<date>/USUV_QC_dropped_<date>.tsv` (dropped by `--qc`).
4. **Many lab rows have host or date `unknown`**: the sheet has no `Species_Latin` or sample date for
   these rows. Fix the sheet and run `--lab --merge` again.
5. **No `local_lab` rows after `--merge`**: check that `usuv_data/local_lab/internal_genomes.fasta` and
   `internal_metadata.csv` exist (run `--lab` first) and that no `usuv_data/lab_*_curated.*` files exist.
6. **Logs**: `usuv_data/logs/` and `usuv_data/usuv_pipeline_<date>.log`.

## Credits
The lab data preparation is based on the R scripts of **Emmanuelle Munger**
(`scripts/emma/240706_SelectSequenceDataNL.R`, `240729_PrepSeqData.R`): the study selection, the
ENA-style metadata fields, the coordinate masking and the N / length rule used in `--qc`.

## Contributing
Contributions are welcome! Please feel free to submit issues or pull requests.

## License
Please refer to the LICENSE file for licensing information.
