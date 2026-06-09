#!/usr/bin/env Rscript
# =============================================================================
# bcx2_liftover_addRSID.R
#
# Lift over BCX2 trans-ethnic meta-analysis summary statistics from hg19 to
# hg38, annotate variants with rsIDs from dbSNP, and write per-trait files
# formatted for LDAK --sum-cors.
#
# Input BCX2 files (hg19, .out.gz) have these key columns:
#   rs_number        : variant ID in "CHR:POS_REF_ALT" format (hg19)
#   reference_allele : reference (A1) allele
#   other_allele     : alternative (A2) allele
#   eaf              : effect allele frequency
#   beta             : effect size estimate
#   se               : standard error
#   p-value          : association p-value
#   n_samples        : sample size
#
# Output per-trait file columns (space-delimited, for LDAK --sum-cors):
#   Predictor  A1  A2  n  BETA  SE  P  rsID
#
# Where Predictor = "CHR:POS" in hg38 coordinates, and rsID is the dbSNP
# identifier (NA if not found). LDAK uses Predictor for allele matching with
# tagging files, while rsID is provided for cross-referencing.
#
# Usage:
#   Rscript bcx2_liftover_addRSID.R \
#     --input-dir  /path/to/bcx2/ta/ \
#     --output-dir /path/to/output/ \
#     --chain      /path/to/hg19ToHg38.over.chain.gz \
#     --dbsnp      /path/to/dbSNP/common_all_20180418_first5cols.txt
#
# Dependencies: data.table, rtracklayer (Bioconductor), optparse
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(rtracklayer)
  library(optparse)
})

# ---------------------------------------------------------------------------
# 1. Parse command-line arguments
# ---------------------------------------------------------------------------
option_list <- list(
  make_option("--input-dir",  type = "character", default = NULL,
              help = "Directory containing BCX2 .out.gz files [required]"),
  make_option("--output-dir", type = "character", default = NULL,
              help = "Directory to write per-trait LDAK-formatted files [required]"),
  make_option("--chain",      type = "character", default = NULL,
              help = "Path to hg19-to-hg38 liftover chain file (may be .gz) [required]"),
  make_option("--dbsnp",      type = "character", default = NULL,
              help = paste("Path to dbSNP common variants file with columns:",
                           "#CHROM POS ID REF ALT [required]")),
  make_option("--pattern",    type = "character", default = "BCX2_.*\\.out\\.gz$",
              help = "Regex pattern to select input files [default: BCX2_.*\\.out\\.gz$]"),
  make_option("--suffix",     type = "character", default = ".hg38.forLDAK.txt",
              help = "Output file suffix appended to trait name [default: .hg38.forLDAK.txt]")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt[["input-dir"]]) || is.null(opt[["output-dir"]]) ||
    is.null(opt[["chain"]])     || is.null(opt[["dbsnp"]])) {
  stop("All four arguments (--input-dir, --output-dir, --chain, --dbsnp) are required.")
}

input_dir  <- normalizePath(opt[["input-dir"]],  mustWork = TRUE)
output_dir <- opt[["output-dir"]]
chain_path <- normalizePath(opt[["chain"]],       mustWork = TRUE)
dbsnp_path <- normalizePath(opt[["dbsnp"]],       mustWork = TRUE)

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# 2. Load liftover chain and dbSNP reference
# ---------------------------------------------------------------------------
message("Loading liftover chain from: ", chain_path)
chain <- import.chain(chain_path)

message("Loading dbSNP reference from: ", dbsnp_path)
# Expected columns: #CHROM  POS  ID  REF  ALT
dbsnp <- fread(dbsnp_path, sep = "\t",
               col.names = c("CHROM", "POS", "rsID", "REF", "ALT"),
               select = 1:5)

# Build a lookup key "CHR:POS" for fast joining
dbsnp[, lookup_key := paste0(CHROM, ":", POS)]

# Keep only the first rsID per position (some positions have multiple entries)
dbsnp <- dbsnp[!duplicated(lookup_key)]
message("  Loaded ", nrow(dbsnp), " dbSNP entries.")

# ---------------------------------------------------------------------------
# 3. Identify input files
# ---------------------------------------------------------------------------
input_files <- list.files(input_dir, pattern = opt[["pattern"]], full.names = TRUE)
if (length(input_files) == 0) {
  stop("No files matching pattern '", opt[["pattern"]], "' found in: ", input_dir)
}
message("Found ", length(input_files), " BCX2 trait file(s) to process.")

# ---------------------------------------------------------------------------
# 4. Per-trait processing loop
# ---------------------------------------------------------------------------
for (fpath in input_files) {

  # Derive trait name from filename (e.g. BCX2_BAS_Trans_GWAMA.out.gz -> BAS)
  fname <- basename(fpath)
  trait <- sub("^BCX2_(.+?)_Trans_GWAMA.*", "\\1", fname)
  message("\n--- Processing trait: ", trait, " ---")

  # ---- 4a. Load and select columns ----------------------------------------
  message("  Reading: ", fname)
  ss <- fread(fpath)

  # Confirm required columns are present
  needed <- c("rs_number", "reference_allele", "other_allele",
              "beta", "se", "p-value", "n_samples")
  missing_cols <- setdiff(needed, names(ss))
  if (length(missing_cols) > 0) {
    warning("  Skipping ", trait, ": missing columns: ",
            paste(missing_cols, collapse = ", "))
    next
  }

  ss <- ss[, .(rs_number, A1 = reference_allele, A2 = other_allele,
               BETA = beta, SE = se, P = `p-value`, n = n_samples)]

  # ---- 4b. Parse hg19 chromosome and position from the variant ID ----------
  # rs_number format: "CHR:POS_REF_ALT"
  message("  Parsing hg19 coordinates from variant ID ...")
  parts <- strsplit(as.character(ss$rs_number), "[:_]")

  # Guard: skip rows where parsing fails (< 2 tokens)
  n_tokens <- lengths(parts)
  if (any(n_tokens < 2)) {
    message("  Dropping ", sum(n_tokens < 2),
            " rows with unparseable variant IDs.")
    ss  <- ss[n_tokens >= 2]
    parts <- parts[n_tokens >= 2]
  }

  ss[, chr_hg19 := sapply(parts, `[[`, 1)]
  ss[, pos_hg19 := as.integer(sapply(parts, `[[`, 2))]

  # Remove rows with missing coordinates
  ss <- ss[!is.na(pos_hg19)]

  # ---- 4c. Liftover hg19 -> hg38 ------------------------------------------
  message("  Lifting over ", nrow(ss), " variants to hg38 ...")

  # Build a GRanges object; seqnames must be "chrN" format for the chain
  gr_hg19 <- GRanges(
    seqnames = paste0("chr", ss$chr_hg19),
    ranges   = IRanges(start = ss$pos_hg19, end = ss$pos_hg19),
    idx      = seq_len(nrow(ss))   # carry original row index through liftover
  )

  lifted <- liftOver(gr_hg19, chain)

  # liftOver returns a GRangesList; entries with no mapping are zero-length
  mapped    <- lengths(lifted) == 1   # exactly one hg38 location
  multi_map <- lengths(lifted) > 1    # ambiguous multi-mapping — discard
  unmapped  <- lengths(lifted) == 0

  message("  Mapped: ", sum(mapped), "  |  Multi-mapped (dropped): ",
          sum(multi_map), "  |  Unmapped (dropped): ", sum(unmapped))

  # Extract hg38 positions for uniquely mapped variants
  lifted_gr <- unlist(lifted[mapped])
  hg38_chr  <- sub("^chr", "", as.character(seqnames(lifted_gr)))
  hg38_pos  <- start(lifted_gr)
  orig_idx  <- lifted_gr$idx

  # Subset the summary-stats table to mapped rows and add hg38 coordinates
  ss_hg38 <- ss[orig_idx]
  ss_hg38[, chr_hg38 := hg38_chr]
  ss_hg38[, pos_hg38 := hg38_pos]

  # LDAK Predictor: "CHR:POS" in hg38 (matching tagging file format)
  ss_hg38[, Predictor := paste0(chr_hg38, ":", pos_hg38)]

  # ---- 4d. Annotate rsIDs from dbSNP (hg38) --------------------------------
  message("  Annotating rsIDs from dbSNP ...")
  ss_hg38[, lookup_key := Predictor]
  ss_hg38 <- merge(ss_hg38, dbsnp[, .(lookup_key, rsID)],
                   by = "lookup_key", all.x = TRUE)

  n_annotated <- sum(!is.na(ss_hg38$rsID))
  message("  rsIDs annotated: ", n_annotated, " / ", nrow(ss_hg38),
          " (", round(100 * n_annotated / nrow(ss_hg38), 1), "%)")

  # ---- 4e. Final column selection and output --------------------------------
  out <- ss_hg38[, .(Predictor, A1, A2, n, BETA, SE, P, rsID)]

  # Remove rows with any NA in the LDAK-essential columns
  essential <- c("Predictor", "A1", "A2", "n", "BETA", "SE", "P")
  out <- out[complete.cases(out[, ..essential])]

  out_path <- file.path(output_dir,
                        paste0("BCX2_", trait, "_Trans_GWAMA", opt[["suffix"]]))
  fwrite(out, out_path, sep = " ", na = "NA", quote = FALSE)

  message("  Written ", nrow(out), " variants -> ", basename(out_path))
}

message("\nDone. Output files are in: ", output_dir)
