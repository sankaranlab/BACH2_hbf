#!/usr/bin/env Rscript
################################################################################
# Prepare METAL meta-analysis output for GCTA-COJO analysis
#
# This script converts METAL output to GCTA-COJO input format by:
# 1. Computing beta and SE from Z-scores and sample sizes
# 2. Orienting alleles to match reference
# 3. Filtering by MAF (>=0.01)
#
# Input: METAL meta-analysis output with rsIDs
# Output: GCTA-COJO ready file with columns: SNP, A1, A2, freq, b, se, p, N
#
# NOTE: INPUT_FILE/OUTPUT_FILE below are environment-specific; adapt the
# /path/to/... placeholders to your own system before running.
#
# References:
# - Beta/SE formulas: https://www.biostars.org/p/319584/
# - METAL format: https://www.nature.com/articles/ng.3538
# - LDSC readme: https://ctg.cncr.nl/documents/p1651/readme.txt
################################################################################

library(data.table)
library(tidyverse)
library(janitor)

# Configuration
INPUT_FILE <- "/path/to/METAL_HBF_rsids/together.tsv"
OUTPUT_FILE <- "/path/to/METAL_HBF_rsids/METAL_HBF_all1.chrALL.MAF01.cojoReady"
MAF_THRESHOLD <- 0.01

# Load METAL output
cat("Loading METAL meta-analysis results...\n")
dt <- fread(INPUT_FILE)
cat(sprintf("Loaded %d variants\n", nrow(dt)))

# Extract original reference allele from rs_number (format: chr:pos:ref:alt)
cat("Extracting reference allele from variant ID...\n")
dt <- dt %>% separate(rs_number, into=c(NA, NA, "originalRef", NA), remove=FALSE)

# Filter by allele frequency and standardize alleles
cat("Filtering by allele frequency...\n")
dt2 <- dt %>% 
    filter(Freq1 >= 0.01 & Freq1 <= 0.99) %>%
    mutate(Allele1 = toupper(Allele1),
           Allele2 = toupper(Allele2))
cat(sprintf("After frequency filter: %d variants remain\n", nrow(dt2)))

# Orient alleles to match original reference
# If Allele1 matches reference, use Freq1; otherwise flip frequency and use Allele1 as alt
cat("Orienting alleles to reference...\n")
dt2b <- as.data.frame(t(apply(dt2, 1, function(x) {
    if (x[['originalRef']] == x[['Allele1']]) {
        x[['freq']] <- x[['Freq1']]
        x[['newalt']] <- x[['Allele2']]
    } else {
        x[['freq']] <- 1 - as.numeric(x[['Freq1']])
        x[['newalt']] <- x[['Allele1']]
    }
    x
})))

# Calculate beta and standard error from Z-scores
# Formulas:
#   Beta = Zscore / sqrt(2 * freq * (1 - freq) * (N + Zscore^2))
#   SE = 1 / sqrt(2 * freq * (1 - freq) * (N + Zscore^2))
cat("Computing beta and SE from Z-scores...\n")
dt3 <- dt2b %>% 
    mutate(freq = as.numeric(freq), 
           Zscore = as.numeric(Zscore), 
           Weight = as.numeric(Weight)) %>%
    mutate(b = Zscore / sqrt(2 * freq * (1 - freq) * (Weight + Zscore^2)),
           se = 1 / sqrt(2 * freq * (1 - freq) * (Weight + Zscore^2))) %>%
    rename(SNP = rsID, 
           A1 = newalt, 
           A2 = originalRef, 
           N = Weight, 
           p = "P-value") %>%
    select(SNP, A1, A2, freq, b, se, p, N)

# Final MAF filter
dt4 <- dt3 %>% filter(freq >= MAF_THRESHOLD)
cat(sprintf("Final variant count after MAF >= %.2f filter: %d\n", MAF_THRESHOLD, nrow(dt4)))

# Write output
cat(sprintf("Writing COJO-ready file to: %s\n", OUTPUT_FILE))
fwrite(dt4, OUTPUT_FILE, sep="\t")

cat("\nSummary statistics:\n")
cat(sprintf("  Mean MAF: %.4f\n", mean(dt4$freq)))
cat(sprintf("  Mean N: %.0f\n", mean(dt4$N)))
cat(sprintf("  Variants with p < 5e-8: %d\n", sum(dt4$p < 5e-8)))

cat("\nDone!\n")
