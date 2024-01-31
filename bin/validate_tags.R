#!/usr/bin/env Rscript

## Script to validate tags (barcodes) used during sample multiplexing
## - Tags should be unique
## - Tag names should be unique
## - Tag names must be alphanumeric and must not contain whitespace, dot, comma, semicolon, or dash
## - Sequencing run ID could be present in tag names (before double underscore)
## - Checks the presence of positive and negative controls
## - Estimates number of unqiue tags and their length
## - For dual assymetric tags, 
##     unique barcodes are converted into a "long" format,
##     a biosample tables (`biosamples_asym.csv` and `biosamples_sym.csv`),
##     and file naming scheme (`file_renaming.tsv`) are exported as well

## Usage:
# validate_tags.R \
#   --tags   tags.fasta \
#   --output tags_validated.fasta



cat("Parsing input options and arguments...\n")

suppressPackageStartupMessages(require(optparse))

## Parse arguments
option_list <- list(
  make_option("--tags",   action="store", default=NA,  type='character', help="FASTA file with tags"),
  make_option("--output", action="store", default=NA,  type='character', help="FASTA file with validated tags")
)
opt <- parse_args(OptionParser(option_list=option_list))


## Validation of the required argiments
if(is.na(opt$tags)){
  cat("Input file is not specified!\n", file=stderr())
  stop()
}
if(is.na(opt$output)){
  cat("Output file is not specified!\n", file=stderr())
  stop()
}

## Assign variables
TAGS <- opt$tags
OUTP <- opt$output

## Log assigned variables
cat(paste("FASTA file with tags: ", TAGS, "\n", sep=""))
cat(paste("Output file with validated tags: ", OUTP, "\n", sep=""))

cat("\n")



############## 

cat("Loading R packages...\n")

load_pckg <- function(pkg = "data.table"){
  suppressPackageStartupMessages( library(package = pkg, character.only = TRUE) )
  cat(paste(pkg, packageVersion(pkg), "\n"))
}

load_pckg("data.table")
load_pckg("Biostrings")

cat("\n")


## Load FASTA sequences
cat("..Loading sequence tables\n")
TAGS <- try( readDNAStringSet(filepath = TAGS, format="fasta") )

if("try-error" %in% class(TAGS)){
  cat("Error in reading FASTA file!\n", file=stderr())
  stop(TAGS)
}

cat("Number of records in the file: ", length(TAGS), "\n")


########################
######################## Validate sample names
########################

cat("\n\n===== Validating sample names =====\n\n")

cat("Positive control: ", 
    ifelse(test = any(grepl(pattern = "PosC", x = names(TAGS))), 
           yes = "Present", no = "Absent"), "\n")

cat("Negative control: ", 
    ifelse(test = any(grepl(pattern = "NegC", x = names(TAGS))), 
           yes = "Present", no = "Absent"), "\n")

## Validate names
cat("\nValidating tag names\n")
newnames <- names(TAGS)

cat("..Replacing leading and trailing spaces and tabs\n")
newnames <- trimws(x = newnames, which = "both")

cat("..Replacing duplicated spaces, dashes, dots, commas, or semicolons\n")
newnames <- gsub(pattern = "\\s+", replacement = " ", x = newnames)
newnames <- gsub(pattern = "\\-+", replacement = "-", x = newnames)
newnames <- gsub(pattern = "\\.+", replacement = ".", x = newnames)
newnames <- gsub(pattern = ",+", replacement = ",", x = newnames)
newnames <- gsub(pattern = ";+", replacement = ";", x = newnames)

cat("..Replacing disallowed symbols\n")
newnames <- iconv(newnames, from = "UTF-8", to = "ASCII//TRANSLIT")
newnames <- gsub(pattern = "[^[:alnum:]]", replacement = "_", x = newnames)

cat("..Replacing the second occurrence of double underscore\n")
## need something like `sed 's/__/_/2g'`
newnames <- sub(pattern  = "__", replacement = "TEMPTEMPTEMPTEMP", x = newnames)
newnames <- gsub(pattern = "_+", replacement = "_", x = newnames)
newnames <- sub(pattern  = "TEMPTEMPTEMPTEMP", replacement = "__", x = newnames)

## test:  c("A__B", "B__C__D", "E__F__G__H", "A___3", "B___4___E")


## Find out which names were changed
renamed <- data.table(
    OriginalName = names(TAGS),
    NewName = newnames)

renamed[ , Renamed := OriginalName != NewName ]
renamed <- renamed[ Renamed == TRUE ]

if(nrow(renamed) > 0){
  cat("..The following tag names were corrected:\n")
  print(
    renamed[, .(OriginalName, NewName)],
    nrows = nrow(renamed), trunc.cols = FALSE)
}



## Check tag name uniqness
nuniq <- length(unique(names(TAGS)))
cat("\nAll tag names unique: ", 
    ifelse(test = (length(TAGS) == nuniq), 
           yes = "TRUE", no = "FALSE"), "\n")

if(length(TAGS) != nuniq){
  cat("..Not all tag names are unique!\n")
  cat("..Resolving tag name uniqness by adding sequential numbers to non-unique names\n")

  dups <- unique(names(TAGS)[ which(duplicated(names(TAGS))) ])
  cat("..Number of duplicates: ", length(dups), "\n")
  cat("..Duplicated names: ", paste(dups, collapse = ", "), "\n")

  dtt <- data.table(ID = 1:length(TAGS), TagName = names(TAGS))
  dtt[ , Duplicated := TagName %in% dups ]
  dtt[ , NewName := TagName ]
  dtt[ 
    Duplicated == TRUE, 
    NewName := paste0(TagName, "_", 1:.N),
    by = "TagName" ]
  setorder(x = dtt, ID)

  names(TAGS) <- dtt$NewName

  rm(dtt)
}


## Check run name
TESTRUN <- grepl(pattern = "__", x = names(TAGS))

cat("\nTag names contain sequencing run ID: ", 
    ifelse(test = any(TESTRUN), 
           yes = "TRUE", no = "FALSE"), "\n")

if(any(TESTRUN)){
  
  ## Check that all samples contain RunID
  cat("\nAll samples contain sequencing run ID: ", 
    ifelse(test = sum(TESTRUN) == length(names(TAGS)), 
           yes = "TRUE", no = "FALSE"), "\n")

  ## Check run name uniqness
  dtt <- data.table(TagName = names(TAGS))
  dtt[ , c("RunID", "SampleID") := tstrsplit(x = TagName, split = "__", keep = 1:2) ]

  cat("Number of run IDs in tag names (ideally, should be = 1): ", length(unique(dtt$RunID)), "\n")
  cat("Run IDs detected: ", paste(unique(dtt$RunID), collapse = ", "), "\n")

  rm(dtt)
}


########################
######################## Validate sequences
########################

cat("\n\n===== Validating sequences =====\n\n")


DUAL <- grepl(pattern = "\\.\\.\\.", x = as.character(TAGS))

if(any(DUAL)){
  cat("Barcode type detected: Dual\n")
} else {
  cat("Barcode type detected: Single (or dual symmetric)\n")
}


##### Single tag

if(any(DUAL) == FALSE){

  cat("Tag length: ", paste(unique(width(TAGS)), collapse = ", "), "\n")

  suniq <- length(unique(as.character(TAGS)))
  cat("\nAll tag sequences unique: ", 
    ifelse(test = (length(TAGS) == suniq), 
           yes = "TRUE", no = "FALSE"), "\n")

if(length(TAGS) != suniq){
  cat("..Not all tag sequences are unique!\n")
  cat("..This should be resolved manually!\n")

  dup_name <- unique( names(TAGS)[ duplicated(as.character(TAGS)) ])
  dup_tags <- as.character(TAGS[ dup_name ])

  cat("..Number of duplicated tags: ", length(dup_name), "\n")

  dupss <- TAGS[ TAGS %in% dup_tags ]
  dups <- data.table(
    TagNames = names(dupss),
    Tags = as.character(dupss))

  dup_smr <- dups[ , .(
    TagNames = paste0("[ ", paste(TagNames, collapse = ", "), " ]")
    ),
    by = "Tags"]

    cat("..Duplicates: \n")
    print(dup_smr, nrows = length(TAGS), trunc.cols = FALSE)

    stop("\nPlease fix the tag sequences!\n")
  }

  ## Export FASTA
cat("Exporting validated tags in FASTA format\n")

writeXStringSet(
  x = TAGS,
  filepath = OUTP,
  compress=FALSE,
  format="fasta",
  width=9999)

} # end of single tag



##### Dual tags

if(any(DUAL) == TRUE){

  ## Convert to tabular format
  dtt <- data.table(
    SampleID = names(TAGS),
    Tags     = as.character(TAGS))

  ## Split dual tags
  dtt[ , c("Tag1", "Tag2") := tstrsplit(x = Tags, split = "\\.\\.\\.", keep = c(1,2)) ]

  ## Check if there are any missing tags
  missing_tags <- dtt[ is.na(Tag1) | is.na(Tag2) ]
  if(nrow(missing_tags) > 0){
    cat("WARNING: missing dual tags detected!\n")
    print(missing_tags)
    stop("\nPlease fix the tag sequences!\n")
  }

  cat("..Forward tag length: ", paste(sort(unique(nchar(dtt$Tag1))), collapse = ", "), "\n")
  cat("..Reverse tag length: ", paste(sort(unique(nchar(dtt$Tag2))), collapse = ", "), "\n")

  cat("\n")
  cat("..Number of unique forward tags: ", length(unique(dtt$Tag1)), "\n")
  cat("..Number of unique reverse tags: ", length(unique(dtt$Tag2)), "\n")

  ## Find unique barcodes
  bu <- data.table(Sequence = unique(c(dtt$Tag1, dtt$Tag2)))

  ## Name unique barcodes
  len <- nchar(nrow(bu))
  bu[ , ID := .I ]
  bu[ , ID := sprintf(paste("%0", len, "d", sep = ""), ID) ]
  bu[ , ID := paste0("bc", ID) ]

  ## Convert to FASTA
  seqs <- DNAStringSet(x = bu$Sequence)
  names(seqs) <- bu$ID



  ## Add bacrode IDs
  dtt <- merge(x = dtt, y = bu, by.x = "Tag1", by.y = "Sequence", all.x = TRUE)
  setnames(x = dtt, old = "ID", new = "ID1")

  dtt <- merge(x = dtt, y = bu, by.x = "Tag2", by.y = "Sequence", all.x = TRUE)
  setnames(x = dtt, old = "ID", new = "ID2")

  dtt[ , Barcodes := paste0(ID1, "--", ID2)]


  dtt[ , TagSymmetry := fifelse(Tag1 == Tag2, yes = "symmetric", no = "asymmetric", na = NA) ]
  cat("Number of symmetric tag combinations: ",  sum(dtt$TagSymmetry %in% "symmetric"),  "\n")
  cat("Number of asymmetric tag combinations: ", sum(dtt$TagSymmetry %in% "asymmetric"), "\n")

  ## Validate barcode combination uniqness
  if(nrow(dtt) != length(unique(dtt$Barcodes))){
    cat("WARNING: non-unique barcode combination detected!\n")

    dups <- dtt$Barcodes[ duplicated(dtt$Barcodes) ]
    print( dtt[ Barcodes %in% dups ] )

    stop("\nPlease fix the tag sequences!\n")
  }


  ## Prepare biosample table for LIMA
  # https://lima.how/faq/biosample.html
  cat("\nExporting biosample tables: 'biosamples_sym.csv' and 'biosamples_asym.csv'\n")
  
  res <- data.table(
    Barcodes     = dtt$Barcodes,
    `Bio Sample` = dtt$SampleID,
    TagSymmetry  = dtt$TagSymmetry)

  fwrite(
    x = res[ TagSymmetry %in% "symmetric", .(Barcodes, `Bio Sample`) ] ,
    file = "biosamples_sym.csv",
    quote = FALSE, sep = ",", col.names = TRUE)

  fwrite(
    x = res[ TagSymmetry %in% "asymmetric", .(Barcodes, `Bio Sample`) ] ,
    file = "biosamples_asym.csv",
    quote = FALSE, sep = ",", col.names = TRUE)


  ## File for sample renaming
  res[ , OldName := paste0("lima.", Barcodes, ".fq.gz") ]
  res[ , NewName := paste0(`Bio Sample`, ".fq.gz") ]

  cat("Exporting file naming scheme: 'file_renaming.tsv'\n")

  fwrite(x = res[ , .(OldName, NewName)],
    file = "file_renaming.tsv", quote = F, sep = "\t", col.names = FALSE)

  ## Export unique barcodes
  cat("Exporting unique tags in FASTA format\n")
  
  writeXStringSet(
    x        = seqs,
    filepath = OUTP,
    compress = FALSE,
    format   = "fasta",
    width    = 9999)

} # end of dual tags


cat("\nValidation finished\n")


##################### Session info

cat("\nAll done.\n")
cat("\n")
cat("Session info:\n")
sessionInfo()
cat("\n")

