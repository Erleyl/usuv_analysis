#!/usr/bin/env Rscript
# =============================================================================
# prepare_lab_fasta.R
# Clean and export the local-lab (and partner) USUV sequences as FASTA + metadata
# for the merge step. Normally run through `bash main.sh --lab`
# (scripts/prepare_lab.sh), which also installs the result for merge_lab.sh.
# =============================================================================
#
# CREDITS
#   The study selection (StudiesInclude / Add / CMV_GenExternal filters), the
#   ENA-style metadata fields (isolate_name, geo_loc_name, sample_title,
#   lat_lon, coordinate masking) and the FASTA export were developed by
#   Emmanuelle Munger in scripts/emma/240706_SelectSequenceDataNL.R (and
#   240729_PrepSeqData.R). This script is a restructured version of her work.
#   Her Included_In study lists are kept and can be switched on (step 5); what
#   was added or changed is listed under "CHANGES" below.
#
# TWO EXPORTS
#   ALL      every row of the sheet that contains sequence data, nothing is
#            filtered out. Each row carries flag columns (Emmanuelle's filter,
#            already on NCBI, repeat, ...) so you can see and choose.
#   CURATION the subset sent to the curation pipeline: everything in ALL except
#            rows that already have a GenBank accession / are already on NCBI,
#            the worse of a sample and its repeat, excluded projects, rows tagged
#            "Excluded" in Included_In,
#            and identical duplicates. Emmanuelle's Included_In study list is
#            NOT applied unless requested with --emma-filter.
#   Both files use the same (unique) SeqIDs.
#   Length and N content are NOT filtered here (n_percent is only reported):
#   scripts/qc_filter.sh (main.sh --qc) does that once, after the merge, for lab
#   and NCBI sequences together.
#
# INPUT
#   --input FILE              USUVSequencingOverview_*.csv (lab overview sheet)
#
# OPTIONS
#   --outdir DIR              output directory (default: usuv_data/lab_sequences)
#   --ncbi-metadata FILE      NCBI metadata (e.g. USUV_Merged_*_metadata.csv).
#                             Used to find lab sequences that are ALREADY on NCBI.
#   --ncbi-fasta FILE         NCBI sequences (e.g. usuv_data/curated/usuv_master.fasta).
#                             Flags rows whose sequence is IDENTICAL to an NCBI record
#                             (review column only, nothing is removed).
#   --external FILE=LABEL     partner sheet(s) with the same columns as the lab sheet,
#                             e.g. usuv_data/local_lab/rivm_seqs.csv=RIVM (comma-separated
#                             for several). Their SeqIDs get the prefix LABEL_ (4940 ->
#                             RIVM_4940); column Source says where a row comes from.
#   --emma-filter MODE        Emmanuelle's Included_In study selection (opt-in):
#                               none (default) | routine | phylodynamics | arbo
#   --emma-extra A,B          Included_In tags added to her list, e.g. 2025
#                             (only used together with --emma-filter routine etc.)
#   --ignore-excluded-tags    also send rows tagged Excluded_* / Exlcuded_* / Not80Perc
#                             in Included_In to CURATION (default: they stay out)
#   --exclude-projects A,B    Project values never sent to curation
#                             (default: Culture,MosquitoSurveillance_CMV_GenExternal)
#   --date-tag YYYYMMDD       tag used in output file names (default: today)
#   --keep-repeats            keep BOTH rows when a sample was sequenced twice
#                             (Remark "Repeat"); default: keep the one with less N
#   --include-published       keep rows that already have a GenBank accession in CURATION
#   --keep-ncbi-duplicates    keep rows whose SeqID is already on NCBI in CURATION
#   --help
#
# OUTPUT (in --outdir, prefix USUV_lab_<date-tag>)
#   _all.fasta, _all_metadata.csv      ALL rows with sequence data (+ flag columns)
#   .fasta, _metadata.csv              CURATION subset  -> merge (installed by prepare_lab.sh)
#   _ENA_submission_metadata.csv       ENA-style table for the CURATION subset
#                                      (own lab rows only: partner data is not ours to submit)
#   _selection_audit.tsv               one line per input row: what happened and why
#
# EXAMPLE
#   Normal use (finds the newest sheet, RIVM and the NCBI master files itself):
#     bash main.sh --lab --merge --qc
#   Direct call:
#     Rscript scripts/prepare_lab_fasta.R \
#       --input usuv_data/local_lab/USUVSequencingOverview_ED261009.csv \
#       --external usuv_data/local_lab/rivm_seqs.csv=RIVM \
#       --ncbi-metadata usuv_data/curated/usuv_master_metadata.csv \
#       --ncbi-fasta usuv_data/curated/usuv_master.fasta
#
# CHANGES compared with the original script
#   * Nothing with sequence data is lost: the ALL export keeps every such row.
#   * Partner sheets (e.g. RIVM) can be added with --external; they go through
#     the same cleaning, ID and curation rules as the lab rows.
#   * The Included_In study list is no longer applied by default, so new samples
#     without a study tag are not lost. When switched on (--emma-filter), tags
#     are matched after normalising spaces/underscores/case (the current sheet
#     writes "Oude_Munnink_et_al_2020", the original list has spaces) and the
#     list can be extended with --emma-extra.
#   * Rows already on NCBI are detected (step 4) and kept out of CURATION.
#   * Sequences are cleaned (gaps/whitespace removed, upper case); placeholder
#     text such as "FAILED" in the Sequence column is not sequence data.
#   * Missing, duplicated or unsafe SeqIDs get a safe, unique ID instead of
#     duplicate FASTA headers.
#   * Collection year comes from DeathOrSampleDate (AnalysisDate only as
#     fall-back); AnalysisDate is blank for many rows.
#   * Latitude/Longitude keep their sign (abs() removed) and masked coordinates
#     are also masked in the metadata passed to the curation step.
#   * sample_title has a fall-back for projects missing from the lookup.
#   * Command-line options instead of setwd() and file names edited by hand.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
})

# ---- 0. Options --------------------------------------------------------------

opt <- list(
  "input"                = NA_character_,
  "outdir"               = "usuv_data/lab_sequences",
  "ncbi-metadata"        = NA_character_,
  "ncbi-fasta"           = NA_character_,
  "external"             = NA_character_,
  "emma-filter"          = "none",
  "emma-extra"           = NA_character_,
  "exclude-projects"     = "Culture,MosquitoSurveillance_CMV_GenExternal",
  "date-tag"             = format(Sys.Date(), "%Y%m%d"),
  "keep-repeats"         = FALSE,
  "include-published"    = FALSE,
  "keep-ncbi-duplicates" = FALSE,
  "ignore-excluded-tags" = FALSE,
  "help"                 = FALSE
)
flags <- c("keep-repeats", "include-published", "keep-ncbi-duplicates", "ignore-excluded-tags", "help")

args <- commandArgs(trailingOnly = TRUE)
# In RStudio: edit the line below instead of using the command line.
if (length(args) == 0 && interactive()) {
  args <- c("--input", "usuv_data/local_lab/USUVSequencingOverview_ED261009.csv")
}

i <- 1
while (i <= length(args)) {
  key <- sub("^--", "", args[i])
  if (!startsWith(args[i], "--") || !key %in% names(opt)) stop("Unknown argument: ", args[i])
  if (key %in% flags) { opt[[key]] <- TRUE; i <- i + 1 }
  else                { opt[[key]] <- args[i + 1]; i <- i + 2 }
}
if (opt[["help"]] || is.na(opt[["input"]])) {
  cat("Usage: Rscript prepare_lab_fasta.R --input overview.csv [--ncbi-metadata merged_metadata.csv] [options]\n",
      "See the header of this script for all options.\n")
  quit(status = if (opt[["help"]]) 0 else 1)
}
if (!file.exists(opt[["input"]])) stop("Input file not found: ", opt[["input"]])
dir.create(opt[["outdir"]], showWarnings = FALSE, recursive = TRUE)

split_list <- function(x) if (is.na(x)) character(0) else str_trim(strsplit(x, ",")[[1]])
normalise  <- function(x) str_to_lower(str_replace_all(str_trim(x), "[\\s_]+", "_"))
out_prefix <- file.path(opt[["outdir"]], paste0("USUV_lab_", opt[["date-tag"]]))

# ---- 1. Read the lab overview ------------------------------------------------
# Everything is read as text so that IDs and dates are never changed by R.
# "UTF-8-BOM": the sheet starts with a byte-order mark (Excel export).

main_label <- "EMC_lab"   # value of the Source column for the lab's own sheet

read_sheet <- function(path) {
  x <- read.csv(path, colClasses = "character", na.strings = c("", "NA"),
                check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
  x[] <- lapply(x, str_trim)
  x
}

df <- read_sheet(opt[["input"]])
df$Source <- main_label
orig_cols <- setdiff(names(df), "Source")

needed <- c("RowNumber", "GenBank", "Project", "SeqID", "UniqueID", "Sequence", "Included_In", "Remark",
            "Country", "province", "Location", "Sample", "Species_Latin",
            "DeathOrSampleDate", "AnalysisDate", "Latitude", "Longitude")
if (any(!needed %in% orig_cols)) stop("Missing columns: ", paste(setdiff(needed, orig_cols), collapse = ", "))
message("Read ", nrow(df), " rows from ", opt[["input"]])

# Re-sort on RowNumber and renumber (a sequence was once added by hand with a
# .5 row number) -- as in the original script.
df <- df %>%
  mutate(RowNumber = suppressWarnings(as.numeric(RowNumber))) %>%
  arrange(RowNumber) %>%
  mutate(Source_RowNumber = RowNumber, RowNumber = row_number())

# 1b. Partner sheets (--external FILE=LABEL). They need the same column names as
#     the lab sheet; at least SeqID and Sequence. Missing columns stay empty,
#     columns that the lab sheet does not have are ignored. RowNumber continues
#     after the lab rows (the original one is kept in Source_RowNumber).
for (spec in split_list(opt[["external"]])) {
  parts <- strsplit(spec, "=", fixed = TRUE)[[1]]
  path  <- parts[1]
  label <- if (length(parts) > 1) parts[2] else toupper(sub("_.*$", "", tools::file_path_sans_ext(basename(path))))
  if (!file.exists(path)) stop("External file not found: ", path)
  if (label == main_label) stop("External label must differ from ", main_label)
  ext <- read_sheet(path)
  if (any(!c("SeqID", "Sequence") %in% names(ext))) stop(path, " needs the columns SeqID and Sequence")
  ignored <- setdiff(names(ext), orig_cols)
  for (col in setdiff(orig_cols, names(ext))) ext[[col]] <- NA_character_
  ext <- ext[, orig_cols]
  ext$Source           <- label
  ext$Source_RowNumber <- suppressWarnings(as.numeric(ext$RowNumber))
  ext$RowNumber        <- max(df$RowNumber) + seq_len(nrow(ext))
  df <- bind_rows(df, ext)
  message("External sheet ", label, ": ", nrow(ext), " rows from ", path,
          if (length(ignored)) paste0(" (ignored columns: ", paste(ignored, collapse = ", "), ")") else "")
}

# ---- 2. Clean the sequences --------------------------------------------------
# Remove alignment gaps ("-"), whitespace and anything that is not a letter;
# upper case. Gap-padded sequences (e.g. old Sanger data) otherwise inflate
# the length and break alignment. A row has sequence data when what is left is
# a nucleotide sequence (IUPAC letters); text such as "FAILED" is not.

df <- df %>%
  mutate(
    Sequence     = toupper(str_replace_all(coalesce(Sequence, ""), "[^A-Za-z]", "")),
    Sequence     = na_if(Sequence, ""),
    seq_len      = nchar(coalesce(Sequence, "")),
    non_iupac    = !is.na(Sequence) & str_detect(Sequence, "[^ACGTRYKMSWBDHVN]"),
    has_sequence = !is.na(Sequence) & !non_iupac,
    n_percent    = round(100 * str_count(coalesce(Sequence, ""), "N") / pmax(seq_len, 1), 1)
  )

# ---- 3. Give every row a safe, unique SeqID ----------------------------------
# Spaces and characters such as : , ( ) ; break IQ-TREE, MAFFT and Newick files.
# The cleaned ID replaces SeqID in FASTA AND metadata so that both match; the
# original is kept in SeqID_original. A row with sequence data but no SeqID is
# NOT dropped: it takes its UniqueID, or LABROW<row number>.
# Repeated IDs get _dup2, _dup3 ... (first occurrence keeps its name).

sanitize <- function(x) na_if(str_replace_all(str_replace_all(x, "\\s+", "_"), "[^A-Za-z0-9._-]", ""), "")

# Partner rows get their label as prefix (4940 -> RIVM_4940): bare numbers could
# collide with other IDs and would not say where the sequence comes from.

df <- df %>%
  mutate(
    SeqID_original = SeqID,
    id_prefix      = if_else(Source != main_label, paste0(Source, "_"), ""),
    id_raw         = if_else(is.na(SeqID), NA_character_, paste0(id_prefix, SeqID)),
    id_clean       = sanitize(id_raw),
    id_unique      = sanitize(if_else(is.na(UniqueID), NA_character_, paste0(id_prefix, UniqueID))),
    SeqID          = coalesce(id_clean, id_unique, sprintf("LABROW%04d", RowNumber)),
    seqid_note     = case_when(
      is.na(id_clean) & !is.na(id_unique) ~ "SeqID_missing_used_UniqueID",
      is.na(id_clean)                      ~ "SeqID_missing_used_RowNumber",
      id_clean != id_raw                   ~ "SeqID_cleaned",
      id_prefix != ""                      ~ paste0("SeqID_prefixed_", Source),
      TRUE                                 ~ NA_character_
    )
  ) %>%
  select(-id_prefix, -id_raw, -id_clean, -id_unique)

idx <- which(df$has_sequence)
# same ID AND same sequence = a true repeat of one entry (flagged, see step 5)
df$dup_identical <- FALSE
df$dup_identical[idx[duplicated(paste(df$SeqID[idx], df$Sequence[idx]))]] <- TRUE
n   <- ave(seq_along(idx), df$SeqID[idx], FUN = seq_along)
ren <- idx[n > 1]
rep_row <- df$Remark[ren] %in% "Repeat"
suffix  <- ifelse(rep_row & n[n > 1] == 2, "_repeat", paste0(ifelse(rep_row, "_repeat", "_dup"), n[n > 1]))
df$seqid_note[ren] <- paste0("renamed_duplicate_of_", df$SeqID[ren])
df$SeqID[ren]      <- paste0(df$SeqID[ren], suffix)

# Same sequence under a DIFFERENT row (e.g. lab sheet vs partner sheet): review
# column only, nothing is removed because of it.
first <- match(df$Sequence, df$Sequence)
df$identical_sequence_as <- NA_character_
hit <- df$has_sequence & first != seq_len(nrow(df))
df$identical_sequence_as[hit] <- df$SeqID[first[hit]]

# ---- 4. Check against NCBI (optional) ----------------------------------------
# Erasmus MC submissions carry the lab SeqID as the LAST "/" part of the NCBI
# `name` (e.g. USUV/NL/2017/MEA155 -> MEA155; older ones are just AS201600045).
# So: lab SeqID == last part of NCBI name  => sequence is already on NCBI.
# ncbi_base = same ID after removing a "_suffix" (e.g. AS201600045_B2 vs
# AS201600045): reported for manual review only, never excluded automatically.

df$ncbi_exact <- NA_character_
df$ncbi_base  <- NA_character_
df$genbank_in_ncbi_download <- NA

if (!is.na(opt[["ncbi-metadata"]])) {
  ncbi <- read.csv(opt[["ncbi-metadata"]], colClasses = "character", na.strings = c("", "NA"),
                   check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
  if ("data_source" %in% names(ncbi)) ncbi <- ncbi[ncbi$data_source %in% "ncbi", ]
  ncbi <- ncbi %>%
    mutate(key = toupper(str_trim(str_extract(name, "[^/]+$")))) %>%
    filter(!is.na(key), key != "UNKNOWN")
  acc_by_key <- vapply(split(ncbi$accession, ncbi$key),
                       function(x) paste(unique(x), collapse = ";"), character(1))

  # Only the lab's own IDs follow the NCBI naming (partner IDs such as bare
  # numbers could match by accident), so partner rows are not compared by ID.
  own      <- df$Source == main_label
  lab_key  <- ifelse(own, toupper(df$SeqID_original), NA_character_)
  lab_base <- ifelse(own, toupper(sub("_.*$", "", df$SeqID_original)), NA_character_)
  df$ncbi_exact <- unname(acc_by_key[lab_key])
  df$ncbi_base  <- ifelse(is.na(df$ncbi_exact) & nchar(lab_base) >= 5, unname(acc_by_key[lab_base]), NA_character_)

  # Sanity check in the other direction: lab sheet says "published" -- is that
  # accession really in the NCBI download?
  ncbi_acc <- sub("\\..*$", "", ncbi$accession)
  df$genbank_in_ncbi_download <- ifelse(is.na(df$GenBank), NA, sub("\\..*$", "", df$GenBank) %in% ncbi_acc)

  message("NCBI check: ", sum(!is.na(df$ncbi_exact)), " lab rows have the same ID as an NCBI record; ",
          sum(!is.na(df$ncbi_base)), " more match only after removing a suffix (review).")
} else {
  message("No --ncbi-metadata given: skipping the NCBI duplicate check.")
}

# 4b. Identical sequence on NCBI (optional, works for every row incl. partner data).
#     Only an exact match counts, so it misses sequences that differ by a single
#     base or by N's; it is a review column, nothing is removed because of it.
df$ncbi_identical_seq <- NA_character_
if (!is.na(opt[["ncbi-fasta"]])) {
  fl   <- readLines(opt[["ncbi-fasta"]])
  isH  <- startsWith(fl, ">")
  ids  <- sub("^>\\s*(\\S+).*$", "\\1", fl[isH])
  nseq <- vapply(split(fl[!isH], cumsum(isH)[!isH]),
                 function(x) toupper(gsub("[^A-Za-z]", "", paste(x, collapse = ""))), character(1))
  m <- match(df$Sequence, nseq)
  df$ncbi_identical_seq <- ifelse(df$has_sequence & !is.na(m), ids[m], NA_character_)
  message("NCBI sequence check: ", sum(!is.na(df$ncbi_identical_seq)), " rows are identical to a record in ",
          opt[["ncbi-fasta"]], " (", length(ids), " records).")
}

# ---- 5. Selection for the curation pipeline ----------------------------------
# Nothing here removes a row from the ALL export; it only decides what also goes
# into the CURATION subset. Rows are tested in this order and the FIRST matching
# reason is recorded.
#
# 5a. Emmanuelle Munger's study selection, copied from her script
#     (StudiesInclude: three versions, the last one "Routine" was the active one;
#      Add: three Bakker et al. bird sequences included by accession).
#     Her filter(Species_English != "Human") for the phylodynamics study is not
#     applied here.

emma_lists <- list(
  phylodynamics = c("Koenraadt_Overwintering", "Oude Munnink et al 2020", "MungerE_Arbosurveillance_2016_2022",
                    "MungerE_Phylodynamics_2016_2023", "Giglia et al in press", "Oude Munnink et al 2019"),
  arbo          = c("Oude Munnink et al 2020", "MungerE_Arbosurveillance_2016_2022",
                    "Giglia et al in press", "Oude Munnink et al 2019", "Koenraadt_Overwintering"),
  routine       = c("Koenraadt_Overwintering", "Oude Munnink et al 2020", "MungerE_Arbosurveillance_2016_2022",
                    "MungerE_Phylodynamics_2016_2023", "Giglia et al in press", "Oude Munnink et al 2019", "2024")
)
emma_add <- c("OP921076", "OP921083", "OP921078")   # birds from Bakker et al.

mode <- opt[["emma-filter"]]
if (!mode %in% c(names(emma_lists), "none")) stop("--emma-filter must be one of: ", paste(c(names(emma_lists), "none"), collapse = ", "))

inc_norm <- normalise(coalesce(df$Included_In, ""))
if (mode == "none") {
  df$emma_filter <- TRUE
} else {
  studies <- normalise(c(emma_lists[[mode]], split_list(opt[["emma-extra"]])))
  df$emma_filter <- inc_norm %in% studies | sub("\\..*$", "", coalesce(df$GenBank, "")) %in% emma_add |
    df$Source != main_label   # her study list only describes the lab's own sheet
}

mark <- function(reason, cond, label) {
  reason[is.na(reason) & !is.na(cond) & cond] <- label
  reason
}
reason <- rep(NA_character_, nrow(df))

# rows without sequence data are in neither export
reason <- mark(reason, is.na(df$Sequence),                                "no_sequence")
reason <- mark(reason, df$non_iupac,                                      "not_a_nucleotide_sequence")  # e.g. the text "FAILED"
# 5a. Emmanuelle's study list (only when --emma-filter is not "none"), the
#     explicit exclusion tags in Included_In, and her CMV_GenExternal exclusion
reason <- mark(reason, !df$emma_filter,                                   "not_in_Emma_filter")
if (!opt[["ignore-excluded-tags"]]) {
  reason <- mark(reason, str_detect(inc_norm, "^(excl|exlc|not80)"),      "Included_In_marked_excluded")  # also catches the "Exlcuded_" typo
}
reason <- mark(reason, df$Project %in% split_list(opt[["exclude-projects"]]), "excluded_project")
# 5b. Extra rules for the curation step (not from the original script)
# A sample sequenced twice (Remark "Repeat"): keep the usable sequence with the
# LEAST N, drop the other. If the first attempt failed, the repeat is simply the
# only usable sequence and is kept. --keep-repeats keeps both.
if (!opt[["keep-repeats"]]) {
  grp      <- paste(df$Source, df$SeqID_original)
  cand     <- df$has_sequence & !is.na(df$SeqID_original)
  is_rep   <- df$Remark %in% "Repeat"
  has_rep  <- ave(as.integer(cand & is_rep), grp, FUN = sum) > 0
  n_usable <- ave(as.integer(cand), grp, FUN = sum)
  quality  <- ifelse(cand, df$n_percent, Inf) - ifelse(is_rep, 1e-6, 0)   # tie: prefer the repeat
  position <- ave(quality, grp, FUN = function(x) rank(x, ties.method = "first"))
  reason   <- mark(reason, cand & has_rep & n_usable > 1 & position > 1, "repeat_or_original_with_more_N")
}
if (!opt[["include-published"]])    reason <- mark(reason, !is.na(df$GenBank),       "has_GenBank_accession_in_lab_sheet")
if (!opt[["keep-ncbi-duplicates"]]) reason <- mark(reason, !is.na(df$ncbi_exact),    "already_on_NCBI")
reason <- mark(reason, df$dup_identical,                                  "duplicate_identical_sequence")

df$in_all_export      <- df$has_sequence
df$in_curation_input  <- is.na(reason)
df$curation_exclusion <- reason

# ---- 6. Dates and coordinates -------------------------------------------------
# Date: sample/death date first, analysis date only as fall-back.
# Week dates (2025-W12) and bare years are reduced to what ENA accepts.
# Coordinates: masked rows are rounded to 2 decimals (~1 km), others to 5.
# DeadCaptiveBirds / HumanBloodBankScreening have no usable coordinates.

date_raw <- str_replace_all(coalesce(df$DeathOrSampleDate, df$AnalysisDate, ""), "_", "-")
df <- df %>%
  mutate(
    date_source     = case_when(!is.na(DeathOrSampleDate) ~ "DeathOrSampleDate",
                                !is.na(AnalysisDate)      ~ "AnalysisDate",
                                TRUE                      ~ NA_character_),
    collection_year = str_extract(date_raw, "^\\d{4}"),
    collection_date = case_when(str_detect(date_raw, "^\\d{4}-\\d{2}(-\\d{2})?$") ~ date_raw,
                                TRUE ~ collection_year)
  )

lat <- suppressWarnings(as.numeric(df$Latitude))
lon <- suppressWarnings(as.numeric(df$Longitude))
masked <- !is.na(df$Remark) & df$Remark == "CoordinatesConfidential_Mask!"
# Impossible values (e.g. 5.19848e+16 = 51.9848 after Excel lost the decimal
# separator) cannot be repaired safely: the coordinates are dropped, the
# sequence is kept, and coordinates_note says what happened.
in_range  <- !is.na(lat) & !is.na(lon) & abs(lat) <= 90 & abs(lon) <= 180
df$coordinates_note <- ifelse(!is.na(lat) & !is.na(lon) & !in_range, "out_of_range_in_sheet_removed", NA_character_)
has_coord <- in_range & !(lat %in% 0 & lon %in% 0) &
  !(df$Project %in% c("DeadCaptiveBirds", "HumanBloodBankScreening"))

lat_txt <- ifelse(masked, sprintf("%.2f", abs(lat)), sprintf("%.5f", abs(lat)))
lon_txt <- ifelse(masked, sprintf("%.2f", abs(lon)), sprintf("%.5f", abs(lon)))
df$lat_lon   <- ifelse(has_coord,
                       paste(lat_txt, ifelse(lat >= 0, "N", "S"), lon_txt, ifelse(lon >= 0, "E", "W")),
                       "not provided")
df$Latitude  <- as.character(ifelse(has_coord, ifelse(masked, round(lat, 2), round(lat, 5)), NA))
df$Longitude <- as.character(ifelse(has_coord, ifelse(masked, round(lon, 2), round(lon, 5)), NA))

# ---- 7. ENA-style fields (from the original script) ---------------------------

join_nonempty <- function(sep, ...) {
  m <- cbind(...)
  apply(m, 1, function(x) paste(x[!is.na(x) & x != ""], collapse = sep))
}
title_map <- c(
  DeadCaptiveBirds           = "Host-derived Usutu virus, bird deceased in captivity, the Netherlands",
  DeadWildBirds              = "Host-derived Usutu virus, dead free-ranging bird, the Netherlands",
  HumanBloodBankScreening    = "Host-derived Usutu virus, Human blood donation, the Netherlands",
  LiveWildBirds              = "Host-derived Usutu virus, live free-ranging birds, the Netherlands",
  LiveWildBirds_Dead         = "Host-derived Usutu virus, dead free-ranging birds, the Netherlands",
  MosquitoesBRS              = "Host-derived Usutu virus, mosquitoes at bird ringing sites, the Netherlands",
  MosquitoSurveillanceWinter = "Host-derived Usutu virus, overwintering mosquitoes, the Netherlands"
)
iso <- c(Netherlands = "NL", Germany = "DE", Belgium = "BE", France = "FR", Spain = "ES",
         Mexico = "MX", Italy = "IT", Austria = "AT", Hungary = "HU", `United Kingdom` = "GB",
         `Czech Republic` = "CZ", Switzerland = "CH")

df <- df %>%
  mutate(
    # geo_loc_name: Country:province_Location (captive birds: Country:province; humans: Country)
    geo_loc_name = case_when(
      Project == "DeadCaptiveBirds"       ~ join_nonempty(":", Country, province),
      Project == "HumanBloodBankScreening" ~ coalesce(Country, ""),
      TRUE                                ~ join_nonempty(":", Country, join_nonempty("_", province, Location))
    ),
    sample_title = coalesce(unname(title_map[Project]),
                            paste0("Host-derived Usutu virus, ", coalesce(Country, "unknown country"))),
    samp_mat_type = Sample %>%
      str_replace_all("_", " ") %>%
      str_replace("MosquitoPool", "mosquito pool") %>%
      str_to_lower() %>%
      str_replace("dbs", "dried blood spot"),
    host = Species_Latin %>%
      str_replace_all("_", " ") %>%
      str_replace_all(fixed("Taeniopygia (genus)"), "Taeniopygia sp."),
    isolate_name = paste("USUV", coalesce(unname(iso[Country]), "XX"),
                         coalesce(collection_year, "unknown"), SeqID, sep = "/")
  )

# ---- 8. Write outputs ---------------------------------------------------------

write_fasta <- function(d, file) {
  wrap <- function(s, w = 70) substring(s, seq(1, nchar(s), w), pmin(seq(1, nchar(s), w) + w - 1, nchar(s)))
  writeLines(unlist(Map(function(id, s) c(paste0(">", id), wrap(s)), d$SeqID, d$Sequence), use.names = FALSE), file)
}
helper_cols <- c("Source", "Source_RowNumber", "SeqID_original", "collection_date", "collection_year",
                 "date_source", "isolate_name", "seqid_note", "coordinates_note", "n_percent",
                 "identical_sequence_as", "ncbi_identical_seq")

all_rows <- df %>% filter(in_all_export)
keep     <- df %>% filter(in_curation_input)

# 8a. ALL rows with sequence data -- FASTA header = SeqID, metadata has flag columns
write_fasta(all_rows, paste0(out_prefix, "_all.fasta"))
all_meta <- all_rows %>%
  select(all_of(orig_cols), all_of(helper_cols), emma_filter, in_curation_input, curation_exclusion,
         ncbi_exact, ncbi_base, genbank_in_ncbi_download) %>%
  rename(ncbi_exact_accession = ncbi_exact, ncbi_base_id_match_REVIEW = ncbi_base)
write.csv(all_meta, paste0(out_prefix, "_all_metadata.csv"), row.names = FALSE, na = "")

# 8b. CURATION subset -- what is sent to the merge step. Lab columns
#     (SeqID = cleaned ID, Latitude/Longitude = masked values, Sequence = gap-free)
#     plus helper columns.
write_fasta(keep, paste0(out_prefix, ".fasta"))
meta <- keep %>%
  select(all_of(orig_cols), all_of(helper_cols), ncbi_base) %>%
  rename(ncbi_base_id_match_REVIEW = ncbi_base)
write.csv(meta, paste0(out_prefix, "_metadata.csv"), row.names = FALSE, na = "")

# 8c. ENA-style submission table for the CURATION subset
#     (original column set; Sequence left out, it is in the FASTA)
ena_cols <- c("GenBank", "sample_title", "isolate_name", "SeqID", "UniqueID", "ID",
              "collection_date", "host", "samp_mat_type", "geo_loc_name", "lat_lon",
              "Latitude", "Longitude", "Path", "Sequencing_Date", "Run", "BC", "Seq_Remark", "Included_In")
ena <- keep %>% filter(Source == main_label) %>%   # partner data is not ours to submit
  select(any_of(ena_cols)) %>% rename(any_of(c(DeadBirdIndID = "ID")))
write.csv(ena, paste0(out_prefix, "_ENA_submission_metadata.csv"), row.names = FALSE, na = "")

# 8d. Audit: what happened to every input row
audit <- df %>%
  select(RowNumber, Source, SeqID_original, SeqID, Project, Included_In, GenBank, Remark, seq_len, n_percent,
         has_sequence, emma_filter, in_all_export, in_curation_input, curation_exclusion,
         ncbi_exact, ncbi_base, ncbi_identical_seq, identical_sequence_as, genbank_in_ncbi_download, date_source) %>%
  rename(ncbi_exact_accession = ncbi_exact, ncbi_base_id_match_REVIEW = ncbi_base)
write.table(audit, paste0(out_prefix, "_selection_audit.tsv"), sep = "\t", row.names = FALSE, na = "", quote = FALSE)

# ---- 9. Summary ---------------------------------------------------------------

message("\n==== Summary ====")
message("Input rows:                                   ", nrow(df))
message("With sequence data (ALL export):              ", nrow(all_rows))
message("Sent to curation (CURATION subset):           ", nrow(keep))
message("Not in the curation subset, by reason:")
print(sort(table(df$curation_exclusion), decreasing = TRUE))
if (length(unique(df$Source)) > 1) {
  message("\nBy source (ALL export / CURATION subset):")
  lv <- unique(df$Source)
  print(rbind(ALL = table(factor(all_rows$Source, levels = lv)), CURATION = table(factor(keep$Source, levels = lv))))
}
outside <- all_rows %>% filter(!emma_filter)
if (nrow(outside) > 0) {
  message("\nNot covered by Emmanuelle's filter (in ALL export only), by Included_In:")
  print(sort(table(coalesce(outside$Included_In, "<blank>")), decreasing = TRUE))
}
if (any(!is.na(all_rows$seqid_note)))   message("\nNOTE: ", sum(!is.na(all_rows$seqid_note)), " SeqIDs were cleaned, generated or renamed (column seqid_note)")
if (any(!is.na(keep$ncbi_base)))        message("REVIEW: ", sum(!is.na(keep$ncbi_base)), " curation rows share a base ID with an NCBI record (column ncbi_base_id_match_REVIEW)")
if (any(!is.na(keep$ncbi_identical_seq)))   message("REVIEW: ", sum(!is.na(keep$ncbi_identical_seq)), " curation rows have a sequence IDENTICAL to an NCBI record (column ncbi_identical_seq)")
if (any(!is.na(keep$identical_sequence_as))) message("REVIEW: ", sum(!is.na(keep$identical_sequence_as)), " curation rows have a sequence identical to another row (column identical_sequence_as)")
if (any(!is.na(all_rows$coordinates_note))) message("WARNING: ", sum(!is.na(all_rows$coordinates_note)), " rows have impossible coordinates in the sheet (e.g. 5.19848e+16, probably a lost decimal separator); coordinates removed, sequences kept (column coordinates_note)")
if (any(keep$n_percent > 10))               message("NOTE: ", sum(keep$n_percent > 10), " curation rows have more than 10% N (column n_percent)")
if (any(df$genbank_in_ncbi_download %in% FALSE)) {
  message("REVIEW: ", sum(df$genbank_in_ncbi_download %in% FALSE),
          " lab rows have a GenBank accession that is NOT in the NCBI download (see audit file)")
}
message("Outputs: ", out_prefix, "{_all.fasta,_all_metadata.csv,.fasta,_metadata.csv,_ENA_submission_metadata.csv,_selection_audit.tsv}")
