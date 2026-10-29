#!/usr/bin/env Rscript

## Script to export rRNA parts (SSU, ITS1, 5.8S, ITS2, LSU, full ITS) and full-length sequences
## for the representative sequences of OTUs

# Input:
#   1. OTU sequences (`OTUs.fa.gz`; IDs = SHA1 hashes of full-length primer-trimmed sequences from Step-1)
#   2. Sequence membership (`UC_Pooled.parquet`)
#   3. Sequence tables in long format (`Seqs.parquet`; sequence of the region used for clustering)
#   4. rRNA parts of Step-1 runs (`<RunID>.parts.parquet`, prepared by `rrna_parts_prep.sh`)
#   5. Sequence qualities of Step-1 runs (`<RunID>.quals.parquet`, prepared by `rrna_parts_prep.sh`)
#   6. (optional) LULU-curated OTU table and LULU merging statistics

# Outputs:
#  - FASTA files with rRNA parts (`SSU.fa.gz`, `ITS1.fa.gz`, `5_8S.fa.gz`, `ITS2.fa.gz`, `LSU.fa.gz`, `ITS_Full.fa.gz`)
#  - FASTA file with full-length primer-trimmed sequences (`OTUs_FullSeq.fa.gz`)
#  - Table with rRNA parts and their origin (`rRNA_parts.parquet`)
#  - Summary (`rRNA_parts_summary.tsv`) and description of the outputs (`README_rRNA_parts.txt`)
#  - The same files with the `LULU_` prefix (for LULU-curated OTUs)

## Selection of rRNA parts:
# For each OTU and rRNA part, the part is taken from the first available sequence ("donor"):
#   0. representative sequence of the OTU
#   1. other members of the OTU (most abundant first, then the highest average Phred score)
#   2. Step-1 sequences with the same sequence of the region used for clustering
#      (e.g., sequences merged during dereplication at the sample level in Step-1)
#   3. members of OTUs merged into the OTU by LULU (and Step-1 sequences with the same region sequence)
# Complete regions are preferred; partial regions (ITSx `--partial`) are used only if no complete region is available

## Usage:
# ./extract_rrna_parts.R \
#    --otus      OTUs.fa.gz \
#    --uc        UC_Pooled.parquet \
#    --seqtab    Seqs.parquet \
#    --parts     parts \
#    --quals     quals \
#    --runs      rRNA_runs.tsv \
#    --luluotus  OTU_table_LULU.txt.gz \
#    --lulustats LULU_merging_statistics.txt.gz \
#    --reptype   "OTUs (VSEARCH clustering at 98% similarity)" \
#    --threads   4 \
#    --memory    8GB

############################################## Parse input parameters

## Check time
start_time <- Sys.time()

cat("Parsing input options and arguments:\n")

suppressPackageStartupMessages(require(optparse))

## Parse arguments
option_list <- list(
  make_option("--otus",      action = "store", default = NA, type = "character", help = "OTU sequences (FASTA)"),
  make_option("--uc",        action = "store", default = NA, type = "character", help = "Sequence membership (UC file in Parquet format)"),
  make_option("--seqtab",    action = "store", default = NA, type = "character", help = "Sequence tables in long format (Parquet)"),
  make_option("--parts",     action = "store", default = NA, type = "character", help = "Directory with rRNA parts of Step-1 runs (*.parquet)"),
  make_option("--quals",     action = "store", default = NA, type = "character", help = "Directory with sequence qualities of Step-1 runs (*.parquet)"),
  make_option("--runs",      action = "store", default = NA, type = "character", help = "Table with Step-1 run info (TSV)"),
  make_option("--luluotus",  action = "store", default = NA, type = "character", help = "LULU-curated OTU table (optional)"),
  make_option("--lulustats", action = "store", default = NA, type = "character", help = "LULU merging statistics (optional)"),
  make_option("--reptype",   action = "store", default = NA, type = "character", help = "Type of representative sequences (for the README)"),
  make_option(c("-t", "--threads"), action = "store", default = 4L, type = "integer", help = "Number of CPU threads for data processing, default 4"),
  make_option(c("-m", "--memory"),  action = "store", default = NA, type = "character", help = "Memory limit for DuckDB (e.g., 8GB)")
)
opt <- parse_args(OptionParser(option_list = option_list))

## Function to convert text "NA"s to NA
to_na <- function(x){
  if(x %in% c("NA", "null", "Null")){ x <- NA }
  return(x)
}

## Replaces "null"s from Nextflow with NA
opt <- lapply(X = opt, FUN = to_na)


## Validation of the required argiments
if(is.na(opt$otus)){
  cat("Input file is not specified: OTU sequences.\n", file = stderr())
  stop()
}
if(is.na(opt$uc)){
  cat("Input file is not specified: UC file is required.\n", file = stderr())
  stop()
}
if(is.na(opt$seqtab)){
  cat("Input file is not specified: sequence tables in Parquet format.\n", file = stderr())
  stop()
}

## Assign variables
OTUS       <- opt$otus
UCF        <- opt$uc
SEQTAB     <- opt$seqtab
PARTS      <- opt$parts
QUALS      <- opt$quals
RUNS       <- opt$runs
LULU_OTUS  <- opt$luluotus
LULU_STATS <- opt$lulustats
REPTYPE    <- opt$reptype
CPUTHREADS <- as.numeric(opt$threads)
MEMORY     <- opt$memory

## Log assigned variables
cat(paste("OTU sequences: ",                  OTUS,       "\n", sep = ""))
cat(paste("UC file (Parquet format): ",       UCF,        "\n", sep = ""))
cat(paste("Sequence tables: ",                SEQTAB,     "\n", sep = ""))
cat(paste("rRNA parts: ",                     PARTS,      "\n", sep = ""))
cat(paste("Sequence qualities: ",             QUALS,      "\n", sep = ""))
cat(paste("Step-1 runs: ",                    RUNS,       "\n", sep = ""))
cat(paste("LULU-curated OTU table: ",         LULU_OTUS,  "\n", sep = ""))
cat(paste("LULU merging statistics: ",        LULU_STATS, "\n", sep = ""))
cat(paste("Representative sequences: ",       REPTYPE,    "\n", sep = ""))
cat(paste("Number of CPU threads to use: ",   CPUTHREADS, "\n", sep = ""))
cat(paste("Memory limit: ",                   MEMORY,     "\n", sep = ""))

cat("\n")


############################################## Data for debugging

# OTUS       <- "OTUs.fa.gz"
# UCF        <- "UC_Pooled.parquet"
# SEQTAB     <- "Seqs.parquet"
# PARTS      <- "parts"
# QUALS      <- "quals"
# RUNS       <- "rRNA_runs.tsv"
# LULU_OTUS  <- "OTU_table_LULU.txt.gz"
# LULU_STATS <- "LULU_merging_statistics.txt.gz"
# REPTYPE    <- NA
# CPUTHREADS <- 4
# MEMORY     <- NA


############################################## Load packages and helpers

cat("Loading R packages:\n")

load_pckg <- function(pkg = "data.table"){
  suppressPackageStartupMessages(library(package = pkg, character.only = TRUE))
  cat(paste(pkg, packageVersion(pkg), "\n"))
}

load_pckg("data.table")
load_pckg("DBI")
load_pckg("duckdb")

cat("\n")

## rRNA parts (column name, output file name, partial-region column)
PART_COLS  <- c("SSU", "ITS1", "S58", "ITS2", "LSU", "ITS")
PART_FILES <- c(SSU = "SSU", ITS1 = "ITS1", S58 = "5_8S", ITS2 = "ITS2", LSU = "LSU", ITS = "ITS_Full")
PART_ALL   <- c(PART_COLS, paste0(setdiff(PART_COLS, "ITS"), "_partial"))

## Donor types
SOURCES <- c("representative", "member", "region_identical", "lulu_member")

## Messages for the README
WARNINGS <- character(0)
add_warning <- function(...){
  msg <- paste0(...)
  cat("WARNING: ", msg, "\n", sep = "")
  WARNINGS <<- c(WARNINGS, msg)
}

## List Parquet files in a directory
list_parquet <- function(x){
  if(is.na(x) || !dir.exists(x)){ return(character(0)) }
  res <- list.files(x, pattern = "\\.parquet$", full.names = TRUE)
  res[ file.size(res) > 0 ]
}

## SQL string literal with a list of files
sql_files <- function(x){
  paste0("[", paste0("'", x, "'", collapse = ", "), "]")
}

## Write FASTA file
write_fasta <- function(ids, seqs, file){
  if(length(ids) == 0){ return(invisible(NULL)) }
  fwrite(
    x = list(c(rbind(paste0(">", ids), seqs))),
    file = file, col.names = FALSE, quote = FALSE, compress = "gzip")
}


############################################## Load the data

cat("..Initializing DuckDB\n")
con <- dbConnect(duckdb(), dbdir = ":memory:")
run_sql <- function(sql){ invisible(dbExecute(con, sql)) }
run_sql(sprintf("SET threads TO %d;", CPUTHREADS))
run_sql("SET preserve_insertion_order = false;")
run_sql("SET temp_directory = 'duckdb_tmp';")
if(!is.na(MEMORY)){
  run_sql(sprintf("SET memory_limit = '%s';", MEMORY))
}

## OTU sequences
cat("..Loading OTU sequences\n")
if(file.size(OTUS) > 0){
  OTU <- fread(
    cmd = sprintf("seqkit fx2tab -i %s | cut -f1,2", shQuote(OTUS)),
    header = FALSE, sep = "\t", col.names = c("OTU", "RepSeq"), colClasses = "character")
} else {
  OTU <- data.table(OTU = character(0), RepSeq = character(0))
}
OTU[ , OTU := tstrsplit(OTU, split = ";", keep = 1) ]
OTU[ , Ord := .I ]
cat("...Number of OTUs: ", nrow(OTU), "\n")
duckdb_register(con, "otus_df", OTU)
run_sql("CREATE TABLE otus AS SELECT * FROM otus_df")

## Sequence membership
cat("..Loading sequence membership\n")
run_sql(sprintf("CREATE VIEW uc_raw AS SELECT * FROM read_parquet('%s')", UCF))
uc_cols <- dbGetQuery(con, "DESCRIBE uc_raw")$column_name
find_col <- function(x){ uc_cols[ tolower(uc_cols) == tolower(x) ][1] }

if(!is.na(find_col("SeqID")) && !is.na(find_col("OTU"))){
  uc_query  <- find_col("SeqID")
  uc_target <- find_col("OTU")
} else if(!is.na(find_col("Query")) && !is.na(find_col("Target"))){
  uc_query  <- find_col("Query")
  uc_target <- find_col("Target")
} else {
  stop("Unexpected columns in the UC file: ", paste(uc_cols, collapse = ", "))
}
cat("...Using columns: ", uc_query, " -> ", uc_target, "\n")

## Members of OTUs (+ representative sequences, which are always members of their OTUs)
run_sql(sprintf("
  CREATE TABLE mem AS
  SELECT SeqID, OTU, CASE WHEN SeqID = OTU THEN 0 ELSE 1 END AS Tier
  FROM (
    SELECT DISTINCT split_part(\"%s\", ';', 1) AS SeqID, split_part(\"%s\", ';', 1) AS OTU
    FROM uc_raw
    WHERE \"%s\" IS NOT NULL
    UNION
    SELECT OTU AS SeqID, OTU FROM otus
  )
  WHERE OTU IN (SELECT OTU FROM otus)",
  uc_query, uc_target, uc_target))

cat("...Number of member sequences: ", dbGetQuery(con, "SELECT count(*) AS n FROM mem")$n, "\n")

## Sequences of the region used for clustering
cat("..Loading sequence tables\n")
run_sql(sprintf("
  CREATE TABLE seqs AS
  SELECT DISTINCT SeqID, upper(Sequence) AS RegionSeq
  FROM read_parquet('%s')
  WHERE SeqID IN (SELECT SeqID FROM mem)",
  SEQTAB))

## Step-1 runs
if(!is.na(RUNS) && file.exists(RUNS) && file.size(RUNS) > 0){
  RUNTAB <- fread(RUNS, sep = "\t", colClasses = "character")
} else {
  RUNTAB <- data.table()
}

## rRNA parts
parts_files <- list_parquet(PARTS)
cat("..Loading rRNA parts: ", length(parts_files), " file(s)\n")

run_sql(paste0(
  "CREATE TABLE cand_raw (RunID VARCHAR, SeqID VARCHAR, Abundance BIGINT, RegionKey VARCHAR, ",
  paste0(PART_ALL, " VARCHAR", collapse = ", "), ")"))

if(length(parts_files) > 0){

  run_sql(sprintf("
    CREATE VIEW parts_all AS
    SELECT * FROM read_parquet(%s, union_by_name = true)",
    sql_files(parts_files)))

  ## Keep only the members of OTUs and Step-1 sequences with the same region sequence
  run_sql(paste0("
    INSERT INTO cand_raw BY NAME
    SELECT RunID, SeqID, Abundance, RegionKey, ", paste(PART_ALL, collapse = ", "), "
    FROM parts_all
    WHERE SeqID IN (SELECT SeqID FROM mem)
       OR RegionKey IN (SELECT RegionSeq FROM seqs)"))
}

cat("...Number of candidate sequences: ",
  dbGetQuery(con, "SELECT count(DISTINCT SeqID) AS n FROM cand_raw")$n, "\n")

## Pool rRNA parts across runs
## (the same hash may be present in multiple runs, keep the most frequent part sequence)
run_sql(paste0("
  CREATE TABLE plong AS
  UNPIVOT (SELECT RunID, SeqID, Abundance, ", paste(PART_ALL, collapse = ", "), " FROM cand_raw)
  ON ", paste(PART_ALL, collapse = ", "), "
  INTO NAME Col VALUE Seq"))

run_sql("
  CREATE TABLE pick AS
  SELECT
    SeqID,
    replace(Col, '_partial', '') AS Part,
    Col LIKE '%_partial' AS Partial,
    Seq
  FROM (SELECT Col, SeqID, Seq, count(*) AS N, sum(Abundance) AS S FROM plong GROUP BY ALL)
  QUALIFY row_number() OVER (PARTITION BY Col, SeqID ORDER BY N DESC, S DESC, length(Seq) DESC, Seq) = 1")

run_sql("
  CREATE TABLE abund AS
  SELECT SeqID, sum(Abundance) AS Abundance
  FROM cand_raw
  GROUP BY SeqID")

run_sql("
  CREATE TABLE regkeys AS
  SELECT DISTINCT SeqID, RegionKey
  FROM cand_raw
  WHERE RegionKey IS NOT NULL")

## Sequence qualities (best read across runs)
quals_files <- list_parquet(QUALS)
cat("..Loading sequence qualities: ", length(quals_files), " file(s)\n")

if(length(quals_files) > 0){
  run_sql(sprintf("
    CREATE TABLE qual AS
    WITH q AS (
      SELECT * FROM read_parquet(%s, union_by_name = true)
      WHERE SeqID IN (SELECT SeqID FROM mem UNION SELECT SeqID FROM cand_raw)
    )
    SELECT
      SeqID,
      CAST(sum(NReads) AS BIGINT) AS NReads,
      arg_max(
        {'RunID': RunID, 'SampleID': SampleID, 'ReadID': ReadID,
         'AvgPhredScore': AvgPhredScore, 'MaxEE': MaxEE, 'MEEP': MEEP, 'FullSeq': FullSeq},
        (coalesce(AvgPhredScore, 0), -coalesce(MaxEE, 0), RunID)
      ) AS Best
    FROM q
    GROUP BY SeqID",
    sql_files(quals_files)))
  run_sql("
    CREATE OR REPLACE TABLE qual AS
    SELECT SeqID, NReads, unnest(Best) FROM qual")
} else {
  run_sql("
    CREATE TABLE qual (SeqID VARCHAR, NReads BIGINT, RunID VARCHAR, SampleID VARCHAR, ReadID VARCHAR,
      AvgPhredScore DOUBLE, MaxEE DOUBLE, MEEP DOUBLE, FullSeq VARCHAR)")
}

if(nrow(OTU) > 0 && length(parts_files) == 0){
  add_warning("No rRNA parts found in the Step-1 results, only full-length sequences are exported")
}
if(nrow(OTU) > 0 && length(quals_files) == 0){
  add_warning("No sequence qualities found in the Step-1 results, full-length sequences are not exported")
}


############################################## Select rRNA parts for OTUs

## `members` = table with SeqID, OTU, and Tier (0 = representative, 1 = OTU member, 3 = member of OTU merged by LULU)
## `otutab`  = table with OTU, RepSeq, and Ord
export_parts <- function(members, otutab, prefix = ""){

  cat("\n..Selecting rRNA parts", ifelse(prefix == "", "", paste0(" (", prefix, ")")), "\n")

  ## Donors (members + Step-1 sequences with the same region sequence)
  run_sql(sprintf("
    CREATE OR REPLACE TABLE donors AS
    WITH memseq AS (
      SELECT m.OTU, s.RegionSeq, min(m.Tier) AS Tier
      FROM %s m JOIN seqs s ON s.SeqID = m.SeqID
      GROUP BY ALL
    ),
    regid AS (
      SELECT ms.OTU, r.SeqID, CASE WHEN ms.Tier = 3 THEN 3 ELSE 2 END AS Tier
      FROM regkeys r JOIN memseq ms ON r.RegionKey = ms.RegionSeq
      WHERE r.SeqID NOT IN (SELECT SeqID FROM %s)
    )
    SELECT OTU, SeqID, min(Tier) AS Tier
    FROM (SELECT OTU, SeqID, Tier FROM %s UNION ALL SELECT * FROM regid)
    GROUP BY OTU, SeqID", members, members, members))

  ## The best donor for each OTU and rRNA part
  run_sql("
    CREATE OR REPLACE TABLE best AS
    SELECT
      d.OTU, p.Part, p.Seq, p.Partial, d.SeqID AS DonorSeqID, d.Tier,
      contains(q.FullSeq, p.Seq) AS InFullSeq
    FROM donors d
    JOIN pick p        ON p.SeqID = d.SeqID
    LEFT JOIN abund a  ON a.SeqID = d.SeqID
    LEFT JOIN qual q   ON q.SeqID = d.SeqID
    QUALIFY row_number() OVER (
      PARTITION BY d.OTU, p.Part
      ORDER BY p.Partial, d.Tier, coalesce(a.Abundance, q.NReads, 0) DESC, q.AvgPhredScore DESC NULLS LAST, d.SeqID
    ) = 1")

  ## Wide table (one row per OTU)
  src_case <- paste0("CASE Tier ", paste0("WHEN ", seq_along(SOURCES) - 1, " THEN '", SOURCES, "'", collapse = " "), " END")
  part_cols <- paste0(sapply(PART_COLS, function(p){
    paste0(
      "max(b.Seq)           FILTER (WHERE b.Part = '", p, "') AS ", p, ", ",
      "max(", src_case, ")  FILTER (WHERE b.Part = '", p, "') AS ", p, "_Source, ",
      "max(b.DonorSeqID)    FILTER (WHERE b.Part = '", p, "') AS ", p, "_DonorSeqID, ",
      "bool_or(b.Partial)   FILTER (WHERE b.Part = '", p, "') AS ", p, "_Partial, ",
      "bool_or(b.InFullSeq) FILTER (WHERE b.Part = '", p, "') AS ", p, "_InFullSeq")
  }), collapse = ",\n      ")

  run_sql(sprintf("
    CREATE OR REPLACE TABLE wide AS
    WITH nm AS (SELECT OTU, count(*) AS NMembers FROM %s GROUP BY OTU),
    bw AS (
      SELECT b.OTU,
      %s
      FROM best b
      GROUP BY b.OTU
    )
    SELECT
      o.OTU, o.Ord, coalesce(nm.NMembers, 1) AS NMembers,
      bw.* EXCLUDE (OTU),
      q.FullSeq, length(q.FullSeq) AS FullSeq_Length,
      q.NReads AS FullSeq_NReads, q.RunID AS FullSeq_RunID, q.SampleID AS FullSeq_SampleID, q.ReadID AS FullSeq_ReadID,
      q.AvgPhredScore AS FullSeq_AvgPhredScore, q.MaxEE AS FullSeq_MaxEE, q.MEEP AS FullSeq_MEEP
    FROM %s o
    LEFT JOIN nm  ON nm.OTU = o.OTU
    LEFT JOIN bw  ON bw.OTU = o.OTU
    LEFT JOIN qual q ON q.SeqID = o.OTU",
    members, part_cols, otutab))

  ## Export the table
  run_sql(sprintf("
    COPY (SELECT * EXCLUDE (Ord) FROM wide ORDER BY Ord)
    TO '%srRNA_parts.parquet' (FORMAT PARQUET, COMPRESSION 'ZSTD', COMPRESSION_LEVEL 8)",
    prefix))

  ## Export FASTA files
  for(p in PART_COLS){
    res <- dbGetQuery(con, sprintf("SELECT OTU, %s AS Seq FROM wide WHERE %s IS NOT NULL ORDER BY Ord", p, p))
    write_fasta(res$OTU, res$Seq, paste0(prefix, PART_FILES[[p]], ".fa.gz"))
  }
  res <- dbGetQuery(con, "SELECT OTU, FullSeq AS Seq FROM wide WHERE FullSeq IS NOT NULL ORDER BY Ord")
  write_fasta(res$OTU, res$Seq, paste0(prefix, "OTUs_FullSeq.fa.gz"))

  ## Summary
  smr <- dbGetQuery(con, paste0("
    SELECT Part, Source, count(*) AS OTUs
    FROM (
      SELECT Part, Source FROM (
        UNPIVOT (SELECT OTU, ", paste0("coalesce(", PART_COLS, "_Source, 'not_found') AS ", PART_COLS, collapse = ", "), " FROM wide)
        ON ", paste(PART_COLS, collapse = ", "), " INTO NAME Part VALUE Source)
      UNION ALL
      SELECT 'FullSeq' AS Part, CASE WHEN FullSeq IS NULL THEN 'not_found' ELSE 'representative' END AS Source FROM wide
    )
    GROUP BY Part, Source"))
  setDT(smr)
  smr[ , Part := factor(Part, levels = c(PART_COLS, "FullSeq")) ]
  smr[ , Source := factor(Source, levels = c(SOURCES, "not_found")) ]
  setorder(smr, Part, Source)
  smr[ , Part := PART_FILES[ as.character(Part) ] ]
  smr[ is.na(Part), Part := "FullSeq" ]
  smr[ , Dataset := ifelse(prefix == "", "OTUs", "LULU") ]
  setcolorder(smr, "Dataset")

  ## Warnings
  nofs <- dbGetQuery(con, "SELECT count(*) AS n FROM wide WHERE FullSeq IS NULL")$n
  if(nofs > 0 && length(quals_files) > 0){
    add_warning(prefix, "Full-length sequence not found for ", nofs, " OTU(s)")
  }
  notin <- dbGetQuery(con, "SELECT count(*) AS n FROM best WHERE InFullSeq = false")$n
  if(notin > 0){
    add_warning(prefix, notin, " rRNA part(s) do not match the full-length sequence of the donor (e.g., reverse-complemented donor)")
  }

  print(smr)
  return(smr)
}

SMR <- export_parts(members = "mem", otutab = "otus", prefix = "")


############################################## LULU-curated OTUs

if(!is.na(LULU_OTUS) && file.exists(LULU_OTUS) && file.size(LULU_OTUS) > 0){

  cat("\n..Processing LULU-curated OTUs\n")

  LULU_IDS <- fread(LULU_OTUS, select = 1, colClasses = "character")[[1]]
  LULU_IDS <- tstrsplit(LULU_IDS, split = ";", keep = 1)[[1]]

  ## Accepted merges (query OTU -> parent OTU)
  ## mumu log has no header, the first column = query OTU, the second = parent OTU, the last = status
  if(!is.na(LULU_STATS) && file.exists(LULU_STATS) && file.size(LULU_STATS) > 0){
    LST <- fread(LULU_STATS, header = FALSE, colClasses = "character")
    LST <- LST[ LST[[ ncol(LST) ]] %in% "accepted", .(Query = V1, Parent = V2) ]
    LST <- unique(LST, by = "Query")
  } else {
    LST <- data.table(Query = character(0), Parent = character(0))
  }
  cat("...Number of merged OTUs: ", nrow(LST), "\n")

  ## Resolve chains of merges (A -> B -> C)
  parent <- setNames(LST$Parent, LST$Query)
  final  <- setNames(LST$Parent, LST$Query)
  for(i in 1:100){
    upd <- !final %in% LULU_IDS & final %in% names(parent)
    if(!any(upd)){ break }
    final[upd] <- parent[ final[upd] ]
  }
  if(any(!final %in% LULU_IDS)){
    add_warning("LULU: ", sum(!final %in% LULU_IDS), " merged OTU(s) could not be assigned to the curated OTUs")
  }
  LMAP <- data.table(OTU = c(LULU_IDS, names(final)), LOTU = c(LULU_IDS, unname(final)))
  LMAP <- unique(LMAP[ LOTU %in% LULU_IDS ], by = "OTU")
  duckdb_register(con, "lulu_map", LMAP)

  ## Members of the curated OTUs
  run_sql("
    CREATE TABLE mem_lulu AS
    SELECT
      m.SeqID, l.LOTU AS OTU,
      CASE WHEN m.SeqID = l.LOTU THEN 0 WHEN m.OTU = l.LOTU THEN m.Tier ELSE 3 END AS Tier
    FROM mem m JOIN lulu_map l ON l.OTU = m.OTU")

  run_sql("
    CREATE TABLE otus_lulu AS
    SELECT * FROM otus WHERE OTU IN (SELECT LOTU FROM lulu_map)")

  nmiss <- length(setdiff(LULU_IDS, OTU$OTU))
  if(nmiss > 0){
    add_warning("LULU: ", nmiss, " curated OTU(s) are not present in the OTU sequences")
  }

  SMR <- rbind(SMR, export_parts(members = "mem_lulu", otutab = "otus_lulu", prefix = "LULU_"))
}


############################################## Summary and README

cat("\n..Writing summary\n")
fwrite(SMR, file = "rRNA_parts_summary.tsv", sep = "\t")

cat("..Writing README\n")
rd <- c(
  "rRNA parts of OTU representative sequences",
  "==========================================",
  "",
  paste0("Representative sequences: ", ifelse(is.na(REPTYPE), "OTUs", REPTYPE)),
  "",
  "OTU IDs are SHA1 hashes of full-length primer-trimmed sequences (Step-1).",
  "Therefore, rRNA parts extracted by ITSx (or ITSx2) in Step-1 and the full-length sequences",
  "could be linked to the representative sequence of each OTU.",
  "",
  "Files:",
  "  SSU.fa.gz, ITS1.fa.gz, 5_8S.fa.gz, ITS2.fa.gz, LSU.fa.gz - rRNA parts",
  "  ITS_Full.fa.gz      - full ITS as extracted by ITSx (`full` region; NB. includes the first nucleotide of LSU)",
  "  OTUs_FullSeq.fa.gz  - full-length primer-trimmed sequence of the representative",
  "  rRNA_parts.parquet  - table with rRNA parts and their origin (one row per OTU)",
  "  rRNA_parts_summary.tsv - number of OTUs with rRNA parts, by origin",
  "  (files with the `LULU_` prefix in the `05.LULU/rRNA_parts` directory - the same for LULU-curated OTUs)",
  "",
  "Origin of rRNA parts (`<Part>_Source` column):",
  "  representative   - the part was extracted from the representative sequence of the OTU",
  "  member           - the part was not detected in the representative sequence,",
  "                     and was taken from the most abundant member of the OTU",
  "                     (ties resolved by the average Phred score)",
  "  region_identical - the part was taken from a Step-1 sequence with the same sequence of the region used for clustering",
  "                     (e.g., merged during dereplication at the sample level in Step-1)",
  "  lulu_member      - the part was taken from a member of an OTU merged into this OTU by LULU",
  "  The donor sequence ID is in the `<Part>_DonorSeqID` column.",
  "  Partial regions (ITSx `--partial`) are used only if a complete region is not available (`<Part>_Partial` column).",
  "  `<Part>_InFullSeq` = the part is a substring of the full-length sequence of the donor.",
  "",
  "Full-length sequence (`FullSeq_*` columns):",
  "  The sequence is identical for all reads with the same hash.",
  "  Quality metrics (AvgPhredScore, MaxEE, MEEP) are reported for the read with the highest average Phred score.",
  "  FullSeq_NReads = number of reads with the same full-length sequence (across all runs).",
  ""
)
if(nrow(RUNTAB) > 0){
  rd <- c(rd, "Step-1 runs:", paste0("  ", c(paste(names(RUNTAB), collapse = "\t"), do.call(paste, c(RUNTAB, sep = "\t")))), "")
}
if(length(WARNINGS) > 0){
  rd <- c(rd, "Warnings:", paste0("  ", unique(WARNINGS)), "")
}
writeLines(rd, "README_rRNA_parts.txt")

## Clean up
dbDisconnect(con, shutdown = TRUE)
if(dir.exists("duckdb_tmp")){ unlink("duckdb_tmp", recursive = TRUE) }

cat("\nAll done.\n")
end_time <- Sys.time()
tmm <- as.numeric(difftime(end_time, start_time, units = "min"))
cat("Elapsed time: ", tmm, " minutes\n")
