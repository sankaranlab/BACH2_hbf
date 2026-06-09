#!/usr/bin/env Rscript
# Extended Data Figure 1A, B, C, D: QQ plots for FEMA and MAMA analyses
# This script generates publication-quality QQ plots with command-line configurability
# Usage: Rscript SupFig1A_QQplot.R <input_file> <output_prefix> <analysis_label>

# Load required libraries
library(tidyverse)      # ggplot2, dplyr
library(data.table)     # fread()
library(ggpubr)         # theme_pubr()

# Parse p-value strings safely and return -log10(P) without numeric underflow.
compute_neglog10p <- function(p_values) {
  p_chr <- tolower(trimws(as.character(p_values)))
  p_chr[p_chr %in% c("", "na", "nan", "inf", "-inf")] <- NA_character_

  out <- rep(NA_real_, length(p_chr))

  sci_match <- stringr::str_match(p_chr, "^([+-]?[0-9]*\\.?[0-9]+)(?:e([+-]?[0-9]+))?$")
  is_valid <- !is.na(sci_match[, 1])

  if (any(is_valid)) {
    mantissa <- suppressWarnings(as.numeric(sci_match[is_valid, 2]))
    exponent <- suppressWarnings(as.numeric(sci_match[is_valid, 3]))
    exponent[is.na(exponent)] <- 0

    valid_num <- !is.na(mantissa) & mantissa > 0 & !is.infinite(mantissa) & !is.na(exponent)
    idx <- which(is_valid)[valid_num]
    out[idx] <- -(log10(mantissa[valid_num]) + exponent[valid_num])
  }

  out
}

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  cat("Usage: Rscript SupFig1A_QQplot.R <input_file> <output_prefix> <analysis_label>\n")
  cat("\nExample:\n")
  cat("  Rscript SupFig1A_QQplot.R CloudTanz_METALOUT_inv_1.tbl ./figures/SupFig1A_1 META\n")
  cat("  Rscript SupFig1A_QQplot.R UA_HBF_MAMA_AFR_HBF.res ./figures/SupFig1A_2 MAMA_AFR\n")
  cat("  Rscript SupFig1A_QQplot.R UA_HBF_MAMA_EUR_HBF.res ./figures/SupFig1A_3 MAMA_EUR\n")
  cat("  Rscript SupFig1A_QQplot.R UA_HBF_MAMA_THAI_HBF.res ./figures/SupFig1A_4 MAMA_EAS\n")
  quit(status = 1)
}

input_file <- args[1]
output_prefix <- args[2]
analysis_label <- args[3]

# Check if input file exists
if (!file.exists(input_file)) {
  cat("Error: Input file does not exist:", input_file, "\n")
  quit(status = 1)
}

cat("Loading data from:", input_file, "\n")
cat("Analysis label:", analysis_label, "\n")

# Read header first to choose p-value column and force it to character on full read.
header <- fread(input_file, nrows = 0)

# Determine which p-value column to use
if ("P-value" %in% colnames(header)) {
  pval_col <- "P-value"
} else if ("P" %in% colnames(header)) {
  pval_col <- "P"
} else {
  cat("Error: Could not find p-value column. Expected 'P-value' or 'P'\n")
  quit(status = 1)
}

cat("P-value column:", pval_col, "\n")

# META GWAS files have 'P-value' column, MAMA files have 'P' column
data <- fread(input_file, colClasses = setNames("character", pval_col))

# Data processing: compute robust observed -log10(P) directly from strings.
data$log10p <- compute_neglog10p(data[[pval_col]])
data <- data[is.finite(log10p) & log10p > 0, ]

# Numeric p-values only for lambda (underflowed values become 0 and are excluded).
data$p_for_lambda <- suppressWarnings(as.numeric(data[[pval_col]]))
lambda_data <- data[is.finite(p_for_lambda) & p_for_lambda > 0 & p_for_lambda <= 1, ]

cat("Total variants:", nrow(data), "\n")
cat("Variants usable for lambda:", nrow(lambda_data), "\n")

# Prepare QQ plot data
ci <- 0.95
n  <- nrow(data)
log10Pe <- expression(paste("Expected -log"[10], plain(P)))
log10Po <- expression(paste("Observed -log"[10], plain(P)))

dt3 <- data %>% 
  arrange(desc(log10p)) %>%
  mutate(observed = log10p) %>%
  mutate(expected = -log10(ppoints(n))) %>%
  mutate(clower = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1))) %>% 
  mutate(cupper = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1)))

# Calculate genomic inflation factor (lambda)
if (nrow(lambda_data) > 0) {
  lambda <- median(qchisq(1 - lambda_data$p_for_lambda, 1), na.rm = TRUE) / qchisq(0.5, 1)
} else {
  lambda <- NA_real_
}
cat("Genomic inflation factor (lambda):", formatC(lambda, format = "f", digits = 3), "\n")

# Create QQ plot
qq_plot <- ggplot(dt3) + 
  geom_point(aes(expected, observed), shape = 16, size = 3) +
  geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
  xlab(log10Pe) +
  ylab(log10Po) + 
  scale_fill_viridis_d() + 
  scale_color_viridis_d() + 
  theme_pubr() + 
  theme(legend.position = "none", text = element_text(size = 16)) +
  ggtitle(paste(analysis_label, "- λ =", formatC(lambda, format = "f", digits = 3)))

# Display plot
# print(qq_plot)

# Save plots
for (ext in c(".png")) { #".pdf", 
  cat("Saving plot to:", paste0(output_prefix, ext), "\n")
  if (ext == ".png") {
    ggsave(filename = paste0(output_prefix, ext), 
           plot = qq_plot, 
           width = 4, height = 4, dpi = 350)
  } else {
    ggsave(filename = paste0(output_prefix, ext), 
           plot = qq_plot, 
           width = 4, height = 4)
  } 
}
cat("Done!\n")
