#!/usr/bin/env Rscript
# Figure 1C, D, E: Manhattan plots from MAMA (Multi-Ancestry Meta Analysis)
# This script generates manhattan plots for ancestry-specific analyses
# Usage: Rscript Figure1CDE_MAMA_manhattan.R <input_file> <output_prefix> <ancestry_label> [figwidth] [figheight] [sig_threshold] [--skip-annotation]

# Load required libraries
library(tidyverse)      # ggplot2, dplyr
library(data.table)     # fread()
library(ggtext)         # element_markdown()
# library(ggbreak)        # scale_y_break()
library(ggrepel)       # geom_text_repel()
library(svglite)        # SVG device for ggsave()


ref_genes_path <- '/lab-share/Hem-Sankaran-e2/Public/ref_genomes/human/hg38_ncbiRefSeq_AllGenes_2026UCSC.tsv.gz'

# Parse p-values from strings to preserve ultra-small values (e.g. 1e-400)
# without requiring external high-precision libraries.
compute_neglog10p <- function(p_values) {
  p_clean <- trimws(tolower(as.character(p_values)))

  parse_one <- function(pv) {
    if (is.na(pv) || pv == "" || pv %in% c("na", "nan", "inf", "-inf")) {
      return(NA_real_)
    }

    # Scientific notation: mantissa e exponent
    m <- regexec("^([+-]?[0-9]*\\.?[0-9]+)(?:e([+-]?[0-9]+))$", pv, perl = TRUE)
    parts <- regmatches(pv, m)[[1]]
    if (length(parts) == 3) {
      mantissa <- suppressWarnings(as.numeric(parts[2]))
      exponent <- suppressWarnings(as.numeric(parts[3]))
      if (is.na(mantissa) || is.na(exponent) || mantissa <= 0) {
        return(NA_real_)
      }
      return(-(log10(mantissa) + exponent))
    }

    # Decimal notation fallback (works for non-underflow values)
    p_num <- suppressWarnings(as.numeric(pv))
    if (is.na(p_num) || p_num <= 0 || is.infinite(p_num)) {
      return(NA_real_)
    }
    -log10(p_num)
  }

  vapply(p_clean, parse_one, numeric(1))
}

# Normalize and standardize summary-stat headers from heterogeneous inputs.
# This handles common variants across FEMA/MAMA/cohort post-QC files.
normalize_input_columns <- function(dt) {
  normalize_key <- function(x) gsub("[^a-z0-9]", "", tolower(x))

  orig_names <- colnames(dt)
  if (is.null(orig_names) || length(orig_names) == 0) {
    return(dt)
  }

  canonical_by_key <- c(
    snp = "SNP", chr = "CHR", bp = "BP", a1 = "A1", a2 = "A2",
    af = "AF", beta = "BETA", se = "SE", z = "Z", p = "P", n = "N"
  )

  alias_to_target <- c(
    snp = "SNP", markername = "SNP", rsid = "SNP", variantid = "SNP",
    chr = "CHR", chromosome = "CHR", chrom = "CHR",
    bp = "BP", position = "BP", genpos = "BP", pos = "BP",
    a1 = "A1", allele2 = "A1", ea = "A1", effectallele = "A1",
    a2 = "A2", allele1 = "A2", nea = "A2", otherallele = "A2",
    af = "AF", freq = "AF", freq1 = "AF", eafreq = "AF", afallele2 = "AF",
    beta = "BETA", effect = "BETA", signedsumstat = "BETA",
    se = "SE", stderr = "SE", sebeta = "SE", betase = "SE",
    z = "Z", zscore = "Z", tstat = "Z",
    p = "P", pvalue = "P", pvart = "P", pbeta = "P", log10p = "P",
    n = "N", weight = "N", totaln = "N", neff = "N", norig = "N"
  )

  keys <- normalize_key(orig_names)
  new_names <- orig_names
  taken_targets <- character(0)
  rename_pairs <- character(0)

  # Pass 1: canonical columns, normalize exact case to expected names.
  for (i in seq_along(orig_names)) {
    key <- keys[i]
    if (!(key %in% names(canonical_by_key))) {
      next
    }
    target <- canonical_by_key[[key]]
    if (target %in% taken_targets) {
      next
    }
    if (orig_names[i] != target) {
      new_names[i] <- target
      rename_pairs <- c(rename_pairs, paste0(orig_names[i], "->", target))
    }
    taken_targets <- c(taken_targets, target)
  }

  # Pass 2: alias columns for remaining targets.
  for (i in seq_along(orig_names)) {
    key <- keys[i]
    if (!(key %in% names(alias_to_target))) {
      next
    }
    target <- alias_to_target[[key]]
    if (target %in% taken_targets) {
      next
    }
    if (orig_names[i] != target) {
      new_names[i] <- target
      rename_pairs <- c(rename_pairs, paste0(orig_names[i], "->", target))
    }
    taken_targets <- c(taken_targets, target)
  }

  if (!identical(orig_names, new_names)) {
    setnames(dt, old = orig_names, new = new_names)
  }
  if (length(rename_pairs) > 0) {
    cat("Header normalization applied:", paste(rename_pairs, collapse = ", "), "\n")
  }

  dt
}

# Decide whether P-like input is raw p-value or already -log10(P).
detect_p_scale <- function(col_name) {
  if (is.null(col_name) || is.na(col_name) || col_name == "") {
    return("raw")
  }
  col_key <- gsub("[^a-z0-9]", "", tolower(col_name))
  if (col_key %in% c("log10p", "minuslog10p", "neglog10p")) {
    return("log10")
  }
  return("raw")
}

# Function to identify significant peaks and annotate with nearest gene
annotate_peaks <- function(gwas_data, gene_annotation, sig_log10_threshold = -log10(5e-8), peak_window = 500000) {
  # Filter for significant variants
  sig_variants <- gwas_data %>%
    filter(log10p > sig_log10_threshold) %>%
    arrange(CHR, BP)
  
  if (nrow(sig_variants) == 0) {
    return(data.frame(
      CHR = integer(), BP = numeric(), P_raw = character(), log10p = numeric(),
      SNP = character(), nearest_gene = character(), peak_label = character()
    ))
  }
  
  # Identify peaks by grouping nearby significant variants
  sig_variants <- sig_variants %>%
    group_by(CHR) %>%
    mutate(
      pos_diff = BP - lag(BP, default = 0),
      new_peak = pos_diff > peak_window | row_number() == 1,
      peak_id = cumsum(new_peak)
    ) %>%
    ungroup()
  
  # Get lead variant for each peak (most significant)
  lead_variants <- sig_variants %>%
    group_by(CHR, peak_id) %>%
    arrange(desc(log10p)) %>%
    slice(1) %>%
    ungroup() %>%
    select(CHR, BP, P_raw, SNP, log10p)
  
  # Annotate with nearest gene
  lead_variants$nearest_gene <- sapply(1:nrow(lead_variants), function(i) {
    chr <- lead_variants$CHR[i]
    pos <- lead_variants$BP[i]
    
    # Filter genes on same chromosome
    chr_genes <- gene_annotation %>%
      filter(chrom == paste0("chr", chr))
    
    if (nrow(chr_genes) == 0) {
      return(NA)
    }
    
    # Calculate distance to gene (0 if variant is within gene)
    chr_genes <- chr_genes %>%
      mutate(
        distance = case_when(
          pos >= txStart & pos <= txEnd ~ 0,
          pos < txStart ~ txStart - pos,
          pos > txEnd ~ pos - txEnd
        )
      ) %>%
      arrange(distance)
    
    # Return nearest gene name
    return(chr_genes$name2[1])
  })
  
  # Format output
  peak_annotation <- lead_variants %>%
    mutate(
      peak_label = paste0("chr", CHR, ":", BP)
    ) %>%
    select(CHR, BP, P_raw, log10p, SNP, nearest_gene, peak_label) %>%
    arrange(desc(log10p))
  
  return(peak_annotation)
}


# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
skip_peak_annotation <- "--skip-annotation" %in% args
args <- args[args != "--skip-annotation"]

if (length(args) < 3) {
  cat("Usage: Rscript Figure1CDE_MAMA_manhattan.R <input_file> <output_prefix> <ancestry_label> [sig_threshold] [figwidth] [figheight] [--skip-annotation]\n")
  cat("\nExample:\n")
  cat("  Rscript Figure1CDE_MAMA_manhattan.R UA_HBF_MAMA_EUR_HBF.res ./figures/Figure1D EUR 6 9\n")
  cat("  Rscript Figure1CDE_MAMA_manhattan.R UA_HBF_MAMA_AFR_HBF.res ./figures/Figure1C AFR 5 4\n")
  cat("  Rscript Figure1CDE_MAMA_manhattan.R UA_HBF_MAMA_EUR_HBF.res ./figures/Figure1D EUR --skip-annotation\n")
  quit(status = 1)
}
print(args)
input_file <- args[1]
output_prefix <- args[2]
ancestry_label <- args[3]
# sig_threshold <- if (length(args) > 3) as.numeric(args[4]) else NULL

# if more than 4 arguments, take the next two as figwidth and figheight (optional)
if (length(args) > 4) {
  figwidth <- as.numeric(args[4])
  figheight <- as.numeric(args[5])
} else {
  figwidth <- 10
  figheight <- 6
}
# add another arg to choose the cutoff for annotating peaks
if (length(args) > 5){
  cutoff_arg <- args[6]
  if (cutoff_arg %in% c("sig", "significant", "5e-8")){
    peak_cutoff <- -log10(5e-8)
  } else if (tolower(cutoff_arg) %in% c("sugg", "suggestive", "sug", "1e-6")){
    peak_cutoff <- -log10(1e-6)
  } else{
    cat("Cannot recognize argument ", cutoff_arg, " for peak annotation cutoff.")
    cat("Falling back to default value...")
    peak_cutoff <- -log10(5e-8)
  }
} else {
  peak_cutoff <- -log10(5e-8)
}
print(peak_cutoff)

# if (!is.null(sig_threshold) && (is.na(sig_threshold) || sig_threshold <= 0)) {
#   cat("Error: sig_threshold must be a positive numeric value.\n")
#   quit(status = 1)
# }
if (is.na(figwidth) || is.na(figheight) || figwidth <= 0 || figheight <= 0) {
  cat("Error: figwidth and figheight must be positive numeric values.\n")
  quit(status = 1)
} else {
  cat("Figure dimensions (width x height):", figwidth, "x", figheight, "inches\n")
}

# Check if input file exists
if (!file.exists(input_file)) {
  cat("Error: Input file does not exist:", input_file, "\n")
  quit(status = 1)
}

cat("Loading data from:", input_file, "\n")

# Load results with robust p-column handling for alternate headers
header_only <- fread(input_file, nrows = 0)
header_names <- names(header_only)
cat("header_names=", paste(header_names, collapse = ", "), "\n")
header_keys <- gsub("[^a-z0-9]", "", tolower(header_names))
cat("header_keys=", paste(header_keys, collapse = ", "), "\n")
p_candidates <- c("p", "pvalue", "log10p", "pvart", "pbeta")
p_idx <- which(header_keys %in% p_candidates)
print(p_idx)
cat("Identified p-value column(s):", if (length(p_idx) > 0) paste(header_names[p_idx], collapse = ", ") else "None", "\n")

if (length(p_idx) > 0) {
  p_col <- header_names[p_idx[1]]
  p_scale <- detect_p_scale(p_col)
  cat("Using p-value column:", p_col, "with detected scale:", p_scale, "\n")
} else {
  p_col <- NA_character_
  p_scale <- "raw"
}

if (p_scale == "raw" && !is.na(p_col)) {
  cat("P values in raw scale. Will read them as characters and compute -log10(P) for plotting.\n")
  data <- fread(input_file, colClasses = setNames("character", p_col))
} else {
  cat("Reading input directly (numeric), then applying standardized column mapping.\n")
  data <- fread(input_file)
}

cat("Before normalizing columns, their names:", paste(colnames(data), collapse = ", "), "\n")
data <- normalize_input_columns(data)
cat("After normalizing:", paste(colnames(data), collapse = ", "), "\n")

required_cols <- c("CHR", "BP", "P")
missing_cols <- setdiff(required_cols, names(data))
if (length(missing_cols) > 0) {
  cat("Error: Missing required columns after header normalization:",
      paste(missing_cols, collapse = ", "), "\n")
  cat("Available columns:", paste(names(data), collapse = ", "), "\n")
  quit(status = 1)
}
if (!("SNP" %in% names(data))) {
  data$SNP <- NA_character_
}

# Data processing
data$BP <- as.numeric(data$BP)
data$P_raw <- data$P
if (p_scale == "log10") {
  data$log10p <- suppressWarnings(as.numeric(data$P_raw))
} else {
  data$log10p <- compute_neglog10p(data$P_raw)
}
print(unique(data$CHR))
data$CHR <- as.integer(data$CHR)
print(unique(data$CHR))

# If chr/pos have NAs, remove NAs
if (any(is.na(data$CHR)) || any(is.na(data$BP))) {
  cat("Warning: Missing values detected in CHR or BP columns. Removing rows with NA in these columns.\n")
  data <- data[!is.na(data$CHR) & !is.na(data$BP), ]
  cat("Total variants after removing NA in CHR/BP:", nrow(data), "\n")
}

cat("Total variants before filtering:", nrow(data), "\n")
data <- data[is.finite(data$log10p) & data$log10p > 2,]
cat("Total variants after filtering (P < 0.01):", nrow(data), "\n")

if (nrow(data) == 0) {
  cat("Error: No variants remain after filtering at P < 0.01.\n")
  quit(status = 1)
}

# Define significance threshold
# if (is.null(sig_threshold)) {
#   sig <- 0.05 / nrow(data)
#   cat("Calculated Bonferroni significance threshold:", formatC(sig, format = "e", digits = 2), "\n")
# } else {
#   sig <- sig_threshold
#   cat("Using provided significance threshold:", formatC(sig, format = "e", digits = 2), "\n")
# }
# sig_log10p <- -log10(sig)
suggestive_log10p <- -log10(1e-6)
genomewide_log10p <- -log10(5e-8)

# Identify and annotate significant peaks (optional)
significant_peaks <- data.frame(
  CHR = integer(), BP = numeric(), P_raw = character(), log10p = numeric(),
  SNP = character(), nearest_gene = character(), peak_label = character()
)

if (skip_peak_annotation) {
  cat("Skipping peak annotation (--skip-annotation enabled).\n")
} else {
  # Read in the gene annotation file
  # Header: #bin\tname\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds\tscore\tname2\tcdsStartStat\tcdsEndStat\texonFrames
  ref_genes <- fread(ref_genes_path)
  # Remove all RNA genes (name ~ ^NR / ^XR)
  cat("Total genes before filtering:", nrow(ref_genes), "\n")
  # ref_genes <- ref_genes[!grepl("^NR", ref_genes$name), ]
  ref_genes <- ref_genes[!grepl("^XR", ref_genes$name), ]
  cat("Total genes after filtering RNA genes:", nrow(ref_genes), "\n")

  significant_peaks <- annotate_peaks(data, ref_genes, 
                                      sig_log10_threshold = peak_cutoff)  
  cat("Identified", nrow(significant_peaks), "significant peaks\n")
  
  if (nrow(significant_peaks) > 0) {
    print(significant_peaks)
  }
  # Hardcode MYB --> HBS1L-MYB for consistency
  significant_peaks$nearest_gene[significant_peaks$nearest_gene == "MYB"] <- "HBS1L-MYB"
  significant_peaks$nearest_gene[significant_peaks$nearest_gene == "HBS1L"] <- "HBS1L-MYB"
  # Hardcode HBB cluster: HBB, HBD, HBG1, HBG2 --> "HBB cluster"
  hbb_genes <- c("HBB", "HBD", "HBG1", "HBG2", "HBE1", "BGLT3", "HBBP1",
                 "OR51V1", "OR51B2", "OR51B4", "OR52Z1P",
                 "OR51B5", "OR51B6", "OR52A1", "OR52A5", "OR52E1", "OR52E2",
                 "OR51M1", "OR51I2", "OR51I1", "OR51Q1", "OR52J3", "OR51L1"
                 )
  significant_peaks$nearest_gene[significant_peaks$nearest_gene %in% hbb_genes] <- "HBB cluster"
  
  if (nrow(significant_peaks) > 0) {
    print(significant_peaks)
  }
}

# Prepare data for plotting: calculate cumulative positions
data_cum <- data |>
  group_by(CHR) |>
  summarise(max_bp = max(BP)) |>
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) |>
  select(CHR, bp_add)

data <- data |>
  inner_join(data_cum, by = "CHR") |>
  mutate(bp_cum = BP + bp_add)

# Calculate axis positions (center of each chromosome)
axis_set <- data |>
  group_by(CHR) |>
  summarize(center = mean(bp_cum))

# Calculate y-axis limits for log-scaling the -log10(P) axis
y_values <- c(data$log10p, suggestive_log10p, genomewide_log10p) #, sig_log10p
ymin <- floor(min(y_values, na.rm = TRUE) * 10) / 10
ymax <- ceiling(max(y_values, na.rm = TRUE)) + 1
if (!is.finite(ymin) || ymin <= 0) {
  ymin <- 1
}
if (!is.finite(ymax) || ymax <= ymin) {
  ymax <- max(ymin + 1, 10)
}

cat("Y-axis limits:", ymin, "to", ymax, "\n")

# Define color scheme based on ancestry
if (ancestry_label == "EUR") {
  colors <- c("#919191", "#39a9db") # blue and gray
} else if (ancestry_label == "AFR") {
  colors <- c("#919191", "#00916e") # green and gray
} else if (ancestry_label == "THAI") {
  colors <- c("#919191", "#e3b505") # yellow/brown and gray
} else {
  colors <- c("#919191", "#7B44E4")  # Default colors
}

# Fall back to generic sans when Arial is unavailable in current R session/device.
plot_font <- "Arial"
if (requireNamespace("systemfonts", quietly = TRUE)) {
  available_families <- unique(systemfonts::system_fonts()$family)
  if (!("Arial" %in% available_families)) {
    plot_font <- "sans"
  }
} else {
  plot_font <- "sans"
}

# Create manhattan plot
manhplot <- ggplot(data, aes(
  x = bp_cum, y = log10p,
  color = as_factor(CHR), size = log10p
)) +
  geom_hline(
    yintercept = suggestive_log10p, color = "#0433ff",
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = genomewide_log10p, color = "#39a9db",
    linetype = "dashed"
  ) +
  geom_point() +
  scale_x_continuous(
    labels = axis_set$CHR,
    breaks = axis_set$center
  ) +
  scale_y_log10(
    expand = c(0, 0),
    limits = c(ymin, ymax),
    breaks = scales::breaks_log(n = 6),
    labels = scales::label_number()
  ) +
  # scale_y_break(c(30, 150)) +
  # Keep data coordinates unchanged while preventing marker clipping at panel edges.
  coord_cartesian(clip = "off") +
  scale_color_manual(values = rep(colors, length.out = length(unique(data$CHR)))) +
  scale_size_continuous(range = c(0.5, 2)) +
  labs(
    x = 'Chromosome',
    y = "-log<sub>10</sub>(<i>P</i>)",
    title = ancestry_label
  ) 

# Annotate significant peaks on the plot
if (nrow(significant_peaks) > 0) {
  # Add cumulative positions to significant peaks for plotting
  significant_peaks <- significant_peaks %>%
    inner_join(data_cum, by = "CHR") %>%
    mutate(bp_cum = BP + bp_add)

  # Remove redundant gene names (keep only the most significant variant per gene)
  significant_peaks <- significant_peaks %>%
    group_by(nearest_gene) %>%
    arrange(desc(log10p)) %>%
    slice(1) %>%
    ungroup()

  print(significant_peaks)

  # Build a repel data frame containing all points so labels avoid
  # collisions with the full Manhattan point cloud, not only labeled peaks.
  repel_data <- data %>%
    mutate(nearest_gene = "") %>%
    left_join(
      significant_peaks %>%
        select(CHR, BP, nearest_gene),
      by = c("CHR", "BP"),
      suffix = c("", ".label")
    ) %>%
    mutate(nearest_gene = ifelse(
      !is.na(nearest_gene.label) & nearest_gene.label != "",
      nearest_gene.label,
      nearest_gene
    )) %>%
    select(-nearest_gene.label)

  manhplot <- manhplot +
    geom_text_repel(
      data = repel_data,
      aes(x = bp_cum, y = log10p, label = nearest_gene),
      size = 4, color = "black",
      family = plot_font, fontface = "italic",
      force = 5,
      box.padding = 0.5,
      point.padding = 0.35,
      segment.color = "grey20",
      segment.size = 0.25,
      arrow = arrow(length = unit(0.01, "npc")),
      ylim = c(log10(genomewide_log10p), NA),
      max.overlaps = Inf,
      # direction = "y",
      inherit.aes = FALSE
    )
}

# Annotate the top-scoring SNP with a horizontal guide to y-axis and value label
top_snp <- data %>%
  arrange(desc(log10p)) %>%
  slice(1)

x_axis_pos <- min(data$bp_cum, na.rm = TRUE)
x_range <- diff(range(data$bp_cum, na.rm = TRUE))
label_x <- x_axis_pos - 0.05 * x_range

top_snp$top_value_label <- formatC(top_snp$log10p, format = "f", digits = 2)

manhplot <- manhplot +
  geom_segment(
    data = top_snp,
    aes(x = x_axis_pos, xend = bp_cum, y = log10p, yend = log10p),
    color = "black",
    linetype = "solid",
    linewidth = 0.3,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = top_snp,
    aes(x = label_x, y = log10p, label = top_value_label),
    hjust = 1,
    vjust = 0.35,
    size = 3.2,
    color = "black",
    inherit.aes = FALSE
  )

# Do aesthetics last:
manhplot <- manhplot +
  theme_minimal(base_family = plot_font) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_markdown(size=11, family = plot_font),
    axis.text.x = element_text(size = 6.5, vjust = 0.5, angle = 45, family = plot_font, color = "black"),
    axis.text.y = element_text(size = 12, hjust = 1, family = plot_font, color = "black"),
    plot.title = element_text(size = 12, hjust = 0.5, family = plot_font),
    panel.border = element_rect(color = "gray60", fill = NA, size = 0.15)
    # plot.margin = margin(t = 10, r = 10, b = 10, l = 80, unit = "pt")
  )

# Save the plot
# Update output paths as needed, "svg", "eps"
for (format in c("png", "pdf")) {
  cat("Saving plots to:\n")
  cat("  ", format, ":", paste0(output_prefix, ".", format), "\n")
  # if it's png, set dpi to 350
  if (format == "png") {
    ggsave(paste0(output_prefix, ".", format), width = figwidth, height = figheight, 
           plot = manhplot, dpi = 350)
  } else {
    pdf_device <- if (capabilities("cairo")) cairo_pdf else "pdf"
    ggsave(paste0(output_prefix, ".", format), width = figwidth, height = figheight, 
         plot = manhplot, device = pdf_device)
  }
}
cat("Done!\n")
