#!/usr/bin/env Rscript
# Extended Data Figure 1E, F, G: beta-beta correlation plots for two sets of sumstats
# MAMA input header: SNP	CHR	BP	A1	A2	FREQ	BETA	SE	Z	P	N_EFF	N_ORIG
# This script generates publication-quality beta-beta correlation plots with command-line configurability
# Usage: Rscript SupFig1_betaCorr.R <input_file1> <input_file2> <output_prefix> <suffix1> <suffix2> [--sig-only]

# Load required libraries
library(tidyverse)      # ggplot2, dplyr
library(data.table)     # fread()
library(ggpubr)         # theme_pubr()

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  cat("Usage: Rscript SupFig1_betaCorr.R <input_file1> <input_file2> <output_prefix> <suffix1> <suffix2> [--sig-only]\n")
  cat("\nExample:\n")
  cat("  Rscript SupFig1_betaCorr.R XC_SNPaligned_AFR_HBF.res XC_SNPaligned_EUR_HBF.res ./figures/SupFig1E EUR AFR\n")
  cat("  Rscript SupFig1_betaCorr.R XC_SNPaligned_AFR_HBF.res XC_SNPaligned_EUR_HBF.res ./figures/SupFig1E EUR AFR --sig-only\n")
  cat("\nOptional argument:\n")
  cat("  --sig-only: If provided, keep only SNPs with P < 5e-8 in both datasets.\n")
  quit(status = 1)
}

input_file1 <- args[1]
input_file2 <- args[2]
output_prefix <- args[3]
suffix1 <- args[4]
suffix2 <- args[5]

# Optional flag: whether to plot only genome-wide significant SNPs
sig_only <- FALSE
if (length(args) >= 6) {
  extra_args <- args[6:length(args)]
  if ("--sig-only" %in% extra_args) {
    sig_only <- TRUE
  }

  unknown_args <- setdiff(extra_args, c("--sig-only"))
  if (length(unknown_args) > 0) {
    cat("Error: Unknown option(s):", paste(unknown_args, collapse = " "), "\n")
    cat("Supported option: --sig-only\n")
    quit(status = 1)
  }
}

# Check if input files exist
if (!file.exists(input_file1)) {
  cat("Error: Input file does not exist:", input_file1, "\n")
  quit(status = 1)
}
if (!file.exists(input_file2)) {
  cat("Error: Input file does not exist:", input_file2, "\n")
  quit(status = 1)
}
analysis_label <- paste0( suffix1, " vs ", suffix2)
cat("Loading data from:", input_file2, "\n")
cat("Analysis label:", analysis_label, "\n")
cat("Significant-only filter (P < 5e-8 in both datasets):", sig_only, "\n")

# Load data based on file structure
cat("Loading data from:", input_file1, "\n")
data1 <- fread(input_file1)
cat("Total variants:", nrow(data1), "\n")

cat("Loading data from:", input_file2, "\n")
data2 <- fread(input_file2)
cat("Total variants:", nrow(data2), "\n")

# Merge data1 and data2 by SNP, CHR, and BP; keep BETA
merged_data <- inner_join(data1 %>% select(SNP, CHR, BP, BETA1 = BETA, P1 = P), 
                          data2 %>% select(SNP, CHR, BP, BETA2 = BETA, P2 = P), 
                          by = c("SNP", "CHR", "BP"))

cat("Overlapping variants:", nrow(merged_data), "\n")

if (sig_only) {
  merged_data <- merged_data %>%
    filter(!is.na(P1), !is.na(P2), P1 < 5e-8, P2 < 5e-8)
  cat("Variants after significance filter:", nrow(merged_data), "\n")
}

if (nrow(merged_data) == 0) {
  cat("Error: No variants available for plotting after applying filters.\n")
  quit(status = 1)
}

# Compute R-squared for BETA1 vs BETA2
r2_val <- summary(lm(BETA2 ~ BETA1, data = merged_data))$r.squared
r2_label <- paste0("R^2 == ", formatC(r2_val, format = "f", digits = 3))

# Create beta-beta correlation plot
beta_corr_plot <- ggplot(merged_data, aes(x = BETA1, y = BETA2)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", col = "red") +
  annotate("text", x = Inf, y = -Inf, label = r2_label, parse = TRUE, hjust = 1.1, vjust = -0.5, size = 4) +
  theme_pubr() +
  xlab(paste0("BETA_", suffix1)) +
  ylab(paste0("BETA_", suffix2)) +
  ggtitle(analysis_label)

# Save plots
for (ext in c(".png")) {  #".pdf", 
  cat("Saving plot to:", paste0(output_prefix, ext), "\n")
  if (ext == ".png") {
    ggsave(filename = paste0(output_prefix, ext), 
           plot = beta_corr_plot, 
           width = 4, height = 4, dpi = 350)
  } else {
    ggsave(filename = paste0(output_prefix, ext), 
           plot = beta_corr_plot, 
           width = 4, height = 4)
  }
}
cat("Done!\n")
