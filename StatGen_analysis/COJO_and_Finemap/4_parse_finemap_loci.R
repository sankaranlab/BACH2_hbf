#!/usr/bin/env Rscript

suppressPackageStartupMessages({
	library(data.table)
})

## ---------------------------------------------------------------------------
## Inputs / outputs
## ---------------------------------------------------------------------------

script_path <- sub("^--file=", "", 
                   grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
script_dir  <- if (length(script_path)) dirname(script_path) else getwd()
root <- system2("git", 
                c("-C", script_dir, "rev-parse", "--show-toplevel"), 
                stdout = TRUE)

## read intervals of significant loci
loci_path     <- file.path(
  root,
  "StatGen_analysis/COJO_and_Finemap/EDFig2c_source_data_5e-8.tsv"
)

## concatenated FINEMAP v1.4.2 results for significant loci
finemap_path  <- file.path(
  root,
  "SCAVENGE_analysis/input_data/hbf-META_finemap_20240712_p5e-8.txt"
)

## gene locations
refgene_path  <- file.path(root, "StatGen_analysis/hg38_ncbiRefSeq_AllGenes_2026UCSC.tsv.gz")

## write the number of credible set SNPs per locus
out_loci_path <- file.path(
  root,
  "StatGen_analysis/COJO_and_Finemap/EDFig2c_source_data_5e-8.tsv"
)

## write the functional element that credible-set SNPs fall into
out_snps_path <- file.path(
  root,
  "StatGen_analysis/COJO_and_Finemap/EDFig2d_source_data_nonzero-cs.tsv"
)

## Feature-model parameters (bp)
## within this distance of the TSS -> "TSS"
TSS_FLANK      <- 500
## further upstream of the TSS (beyond TSS_FLANK) -> "promoter"
PROMOTER_FLANK <- 1000

## ---------------------------------------------------------------------------
## Load locus window definitions
## ---------------------------------------------------------------------------
loci <- fread(loci_path)
setnames(loci, tolower(names(loci)))

## ---------------------------------------------------------------------------
## Load FINEMAP variants and keep only credible-set members (CREDIBLE_SET == 0
## means the variant was not assigned to any credible set)
## ---------------------------------------------------------------------------
finemap <- fread(finemap_path)
finemap[, CREDIBLE_SET := suppressWarnings(as.integer(CREDIBLE_SET))]
credset <- finemap[!is.na(CREDIBLE_SET) & CREDIBLE_SET != 0L]
message("Loaded ", nrow(credset),
        " credible-set variants across ",
        uniqueN(credset$CHR), " chromosomes")

## ---------------------------------------------------------------------------
## Assign each variant to a locus window from EDFig2c_source_data.tsv
## ---------------------------------------------------------------------------
assign_locus <- function(c, b) {
  chr <- win_start <- win_end <- NULL
  hit <- loci[chr == c & b >= win_start & b <= win_end]
  if (nrow(hit) == 0) return(NA_character_)
  paste(hit$window_name, collapse = ";")
}

credset[, locus := mapply(assign_locus, CHR, BP)]

n_unassigned <- sum(is.na(credset$locus))
if (n_unassigned > 0) {
  message("Warning: ", n_unassigned,
          " credible-set variant(s) fall outside any defined locus window")
}

## ---------------------------------------------------------------------------
## Load UCSC ncbiRefSeq gene models, restricted to the chromosomes we need
## ---------------------------------------------------------------------------
chroms_needed <- paste0("chr", unique(credset$CHR))

refgene <- fread(refgene_path)
setnames(refgene, "#bin", "bin")
refgene <- refgene[chrom %in% chroms_needed]

## txStart/cdsStart are 0-based (UCSC convention); convert to 1-based coordinates
refgene[, txStart1  := txStart + 1L]
refgene[, cdsStart1 := cdsStart + 1L]

## strand-aware transcription start site (1-based)
refgene[, tss := ifelse(strand == "+", txStart1, txEnd)]

parse_exons <- function(starts, ends) {
  s <- as.numeric(strsplit(starts, ",")[[1]]) + 1L
  e <- as.numeric(strsplit(ends, ",")[[1]])
  list(start = s, end = e)
}

in_any_exon <- function(pos, starts, ends) {
  ex <- parse_exons(starts, ends)
  any(pos >= ex$start & pos <= ex$end)
}

## ---------------------------------------------------------------------------
## Classify a single variant against overlapping/nearby gene models.
## Priority when several transcripts/genes overlap the same position:
## TSS > 3'UTR > exon > intron > promoter, with ties broken by preferring the
## gene named in the variant's locus window and then by proximity to the TSS.
## ---------------------------------------------------------------------------
classify_variant <- function(chr, pos, locus_name) {
  chrom <- dist_tss <- tss <- in_body <- txStart1 <- txEnd <- NULL
  at_tss <- at_promoter <- strand <- locus_match <- name2 <- NULL
  abs_dist_tss <- rank <- NULL
  chrom_full <- paste0("chr", chr)
  cand <- refgene[chrom == chrom_full]
  if (nrow(cand) == 0) {
    return(list(gene = NA_character_, element = "intergenic"))
  }

  cand <- copy(cand)
  # >0 = downstream of TSS, <0 = upstream
  cand[, dist_tss := pos - tss]
  cand[, in_body     := pos >= txStart1 & pos <= txEnd]
  cand[, at_tss       := abs(dist_tss) <= TSS_FLANK]
  cand[, at_promoter  := !in_body & !at_tss &
         ((strand == "+" & dist_tss < 0 & dist_tss >= -PROMOTER_FLANK) |
          (strand == "-" & dist_tss > 0 & dist_tss <= PROMOTER_FLANK))]

  hits <- cand[in_body | at_tss | at_promoter]
  if (nrow(hits) == 0) {
    return(list(gene = NA_character_, element = "intergenic"))
  }

  hits[, element := "intergenic"]
  for (i in seq_len(nrow(hits))) {
    row <- hits[i]
    if (isTRUE(row$at_tss)) {
      hits[i, element := "TSS"]
    } else if (isTRUE(row$in_body)) {
      is_exon <- in_any_exon(pos, row$exonStarts, row$exonEnds)
      if (!is_exon) {
        hits[i, element := "intron"]
      } else if (row$cdsStart1 > row$cdsEnd) {
        ## non-coding transcript: no CDS, so no UTR distinction
        hits[i, element := "exon"]
      } else if (row$strand == "+" && pos > row$cdsEnd) {
        hits[i, element := "3'UTR"]
      } else if (row$strand == "-" && pos < row$cdsStart1) {
        hits[i, element := "3'UTR"]
      } else {
        hits[i, element := "exon"]
      }
    } else if (isTRUE(row$at_promoter)) {
      hits[i, element := "promoter"]
    }
  }

  priority <- c("TSS" = 1, "3'UTR" = 2,
                "exon" = 3, "intron" = 4, "promoter" = 5)
  hits[, rank := priority[element]]

  ## prefer the gene named by the locus window (e.g. "HBS1L-MYB", "CTC1;KRBA2")
  locus_genes <- if (!is.na(locus_name)) {
    unlist(strsplit(locus_name, "[;-]"))
  } else {
    character(0)
  }
  hits[, locus_match := name2 %in% locus_genes]
  hits[, abs_dist_tss := abs(dist_tss)]

  setorder(hits, -locus_match, rank, abs_dist_tss)
  best <- hits[1]
  list(gene = best$name2, element = best$element)
}

message("Classifying ", nrow(credset),
        " variants against RefSeq gene models...")
ann <- rbindlist(lapply(seq_len(nrow(credset)), function(i) {
  as.data.table(classify_variant(
    credset$CHR[i],
    credset$BP[i],
    credset$locus[i]
  ))
}))

credset <- cbind(credset, ann)

## ---------------------------------------------------------------------------
## Count credible-set variants per locus and update EDFig2c_source_data.tsv
## ---------------------------------------------------------------------------
counts <- credset[!is.na(locus), .N, by = locus]
loci[, finemap_credset := counts$N[match(window_name, counts$locus)]]
loci[is.na(finemap_credset), finemap_credset := 0L]

fwrite(loci, out_loci_path, sep = "\t")
message("Updated locus credible-set counts written to: ", out_loci_path)

## ---------------------------------------------------------------------------
## Write per-variant gene/functional-element annotation table (EDFig2d)
## ---------------------------------------------------------------------------
out_cols <- c("locus", "SNP", "CHR", "BP", "A1", "A2",
              "PIP", "CREDIBLE_SET", "gene", "element")
out <- credset[, ..out_cols]
setnames(out, c("locus", "SNP", "CHR", "BP", "A1", "A2",
                "PIP", "CREDIBLE_SET", "gene", "functional_element"))
setorder(out, locus, CHR, BP)

fwrite(out, out_snps_path, sep = "\t")
message("Wrote variant-level annotations to: ", out_snps_path)
