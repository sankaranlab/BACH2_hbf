#!/usr/bin/env Rscript
# Figure 1B: Manhattan plot from META GWAS with Tanzania cohort
# This script generates the main manhattan plot for the HbF BACH2 publication

# Load required libraries
library(tidyverse)      # ggplot2, dplyr (ggplot, group_by, summarise, mutate, inner_join, filter, select, pull)
library(data.table)     # fread()
library(ggtext)         # element_markdown()
library(ggbreak)        # scale_y_break()
library(ggrepel)       # geom_text_repel()
library(svglite)        # SVG device for ggsave()

# Take in output_prefix and figsize from command line arguments (if any)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >=4) {
  meta_file <- args[1]
  output_prefix <- args[2]
  figwidth <- as.numeric(args[3])
  figheight <- as.numeric(args[4])
} else {
  # give default values
  meta_file <- '/lab-share/Hem-Sankaran-e2/Public/projects/xhcheng/HbF/cohort_sumstats/FEMA/CloudTanz_METALOUT_inv_1.tbl'
  output_prefix <- "Figure1B_manhattan_plot"
  figwidth <- 10
  figheight <- 6
}

# Load META GWAS results
cat(format(Sys.time(), "%Y%m%d_%H%M%S"), " Loading META GWAS results...\n")
# This file is the direct output of the METAL run and is 2.5G in size
# META <- fread('/broad/sankaranlab/arora/hbf/reviews/FEMA/with_tanz_hbf/CloudTanz_METALOUT_inv_1.tbl')
## Local file path on E3:
# META <- fread('/lab-share/Hem-Sankaran-e2/Public/projects/xhcheng/HbF/cohort_sumstats/FEMA/CloudTanz_METALOUT_inv_1.tbl')
META <- fread(meta_file, colClasses = c("P-value" = "character"))
cat(format(Sys.time(), "%Y%m%d_%H%M%S"), " META GWAS results loaded. Total variants:", nrow(META), "\n")

# If "CHR" and "GENPOS" colnames exist, rename to "Chromosome" and "Position" for consistency
if (all(c("CHR", "GENPOS") %in% colnames(META))) {
  META <- META %>%
    rename(Chromosome = CHR, Position = GENPOS)
}
# Strip "chr" prefix if present, and convert to ordered factor for correct sorting
# META$Chromosome <- as.character(META$Chromosome) %>%
#   str_remove("^chr") %>%
#   factor(levels = c(as.character(1:22), "X", "Y", "M"), ordered = TRUE)
META$Chromosome <- as.numeric(as.character(META$Chromosome))
META$Position <- as.numeric(META$Position)

# Data processing
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

META$P_value_raw <- META$`P-value`
META$log10p <- compute_neglog10p(META$P_value_raw)
META <- META[!is.na(META$Chromosome) & !is.na(META$Position), ]
META <- META[!is.na(META$log10p) & is.finite(META$log10p), ]

# Filter to p < 0.001 to avoid cluttering the plot with too many dots
cat("Total variants before p-value filtering:", nrow(META), "\n")
META <- META[META$log10p > 3,]
# META <- META[META$`P-value` < 0.01,]
cat("Total variants after p-value filtering:", nrow(META), "\n")

# Calculate significance threshold
# sig <- 0.05 / nrow(META)
sig <- 5e-8
suggestive_sig <- 1e-6
sig_log10p <- -log10(sig)
suggestive_log10p <- -log10(suggestive_sig)


# Read in the gene annotation file
# Header: #bin	name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	score	name2	cdsStartStat	cdsEndStat	exonFrames
ref_genes <- fread('/lab-share/Hem-Sankaran-e2/Public/ref_genomes/human/hg38_ncbiRefSeq_AllGenes_2026UCSC.tsv.gz')
# Remove all RNA genes (name ~ ^NR)
cat("Total genes before filtering:", nrow(ref_genes), "\n")
ref_genes <- ref_genes[!grepl("^NR", ref_genes$name), ]
ref_genes <- ref_genes[!grepl("^XR", ref_genes$name), ]
cat("Total genes after filtering RNA genes:", nrow(ref_genes), "\n")

# Function to identify significant peaks and annotate with nearest gene
annotate_peaks <- function(gwas_data, gene_annotation, sig_log10_threshold = -log10(5e-8), peak_window = 500000) {
  # Filter for significant variants
  sig_variants <- gwas_data %>%
    filter(log10p > sig_log10_threshold) %>%
    arrange(Chromosome, Position)
  
  if (nrow(sig_variants) == 0) {
    return(data.frame())
  }
  
  # Identify peaks by grouping nearby significant variants
  sig_variants <- sig_variants %>%
    group_by(Chromosome) %>%
    mutate(
      pos_diff = Position - lag(Position, default = 0),
      new_peak = pos_diff > peak_window | row_number() == 1,
      peak_id = cumsum(new_peak)
    ) %>%
    ungroup()
  
  # Get lead variant for each peak (most significant)
  lead_variants <- sig_variants %>%
    group_by(Chromosome, peak_id) %>%
    arrange(desc(log10p)) %>%
    slice(1) %>%
    ungroup() %>%
    select(Chromosome, Position, MarkerName, P_value_raw, log10p)
  
  # Annotate with nearest gene
  lead_variants$nearest_gene <- sapply(1:nrow(lead_variants), function(i) {
    chr <- lead_variants$Chromosome[i]
    pos <- lead_variants$Position[i]
    
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
      peak_label = paste0("chr", Chromosome, ":", Position)
    ) %>%
    select(Chromosome, Position, P_value_raw, log10p, MarkerName, nearest_gene, peak_label) %>%
    arrange(desc(log10p))
  
  return(peak_annotation)
}

# Identify and annotate significant peaks
significant_peaks <- annotate_peaks(META, ref_genes, sig_log10_threshold = sig_log10p)
cat("Identified", nrow(significant_peaks), "significant peaks\n")
# Hardcode MYB --> HBS1L-MYB for consistency
significant_peaks$nearest_gene[significant_peaks$nearest_gene == "MYB"] <- "HBS1L-MYB"
if (nrow(significant_peaks) > 0) {
  print(significant_peaks)
}


# Prepare data for plotting: calculate cumulative positions
data_cum <- META |>
  group_by(Chromosome) |>
  summarise(max_bp = max(Position)) |>
  arrange(Chromosome) |>
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) |>
  select(Chromosome, bp_add)

META <- META |>
  inner_join(data_cum, by = "Chromosome") |>
  mutate(bp_cum = Position + bp_add)

# Calculate axis positions (center of each chromosome)
axis_set <- META |>
  group_by(Chromosome) |>
  summarize(center = mean(bp_cum)) |>
  arrange(Chromosome)  # factors sort by level order

# Calculate y-axis limit based on most significant variant
ymin <- min(c(META$log10p, sig_log10p, suggestive_log10p), na.rm = TRUE)
ymin <- floor(ymin * 10) / 10
ylim <- META |>
  summarize(ylim = ceiling(max(log10p, na.rm = TRUE)) + 8) |>
  pull(ylim)
if (!is.finite(ymin) || ymin <= 0) {
  ymin <- 1
}
if (!is.finite(ylim) || ylim <= ymin) {
  ylim <- max(ymin + 1, 10)
}

# Create manhattan plot
manhplot <- ggplot(META, aes(
  x = bp_cum, y = log10p,
  color = as_factor(Chromosome), size = log10p
)) +
  geom_hline(
    yintercept = sig_log10p, color = "darkred",
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = suggestive_log10p, color = "darkblue",
    linetype = "dashed"
  ) +
  geom_point() +
  scale_x_continuous(
    label = axis_set$Chromosome,
    breaks = axis_set$center
  ) +
  scale_y_log10(
    expand = c(0, 0),
    limits = c(ymin, ylim),
    breaks = scales::breaks_log(n = 6),
    labels = scales::label_number()
  ) +
  scale_color_manual(values = rep(
    c("#6C757D", "#7B44E4"),
    unique(length(axis_set$Chromosome))
  )) +
  scale_size_continuous(range = c(0.5, 2)) +
  # Keep data coordinates unchanged while preventing marker clipping at panel edges.
  coord_cartesian(clip = "off") +
  labs(
    x = 'Chromosome',
    y = "-log<sub>10</sub>(p)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_markdown(),
    axis.text.x = element_text(size = 8, vjust = 0.5), # angle = 60, 
    axis.line.y.left = element_line(color = "black", linewidth = 0.5),
    axis.ticks.y.left = element_line(color = "black", linewidth = 0.5)
  )

if (ymin < 20 && ylim > 180) {
  manhplot <- manhplot + scale_y_break(c(20, 180))
}

# Annotate significant peaks on the plot
if (nrow(significant_peaks) > 0) {
  # Add cumulative positions to significant peaks for plotting
  significant_peaks <- significant_peaks %>%
    inner_join(data_cum, by = "Chromosome") %>%
    mutate(bp_cum = Position + bp_add)
    # Remove redundant gene names (keep only the most significant variant per gene)
    significant_peaks <- significant_peaks %>%
      group_by(nearest_gene) %>%
      arrange(desc(log10p)) %>%
      slice(1) %>%
      ungroup()
    print(significant_peaks)
    
    # Separate peaks by their position relative to the y-axis break
    peaks_lower <- significant_peaks %>% filter(log10p <= 30)
    peaks_upper <- significant_peaks %>% filter(log10p >= 150)
  
  manhplot <- manhplot +
    geom_text(  # _repel
      data = peaks_lower,
      aes(x = bp_cum, y = log10p, label = nearest_gene),
      size = 3.5, color="black", hjust=-0.1, vjust=-0.1,
      # box.padding = 0.5,
      # point.padding = 0.3,
      # segment.color = "grey50",
      # segment.size = 0.3,
      # max.overlaps = 20,
      # ylim = c(NA, 30),
      # inherit.aes = FALSE
    ) +
    geom_text(  #_repel
      data = peaks_upper,
      aes(x = bp_cum, y = log10p, label = nearest_gene),
      size = 3.5, color="black", hjust=-0.1, vjust=1,
      # box.padding = 0.5,
      # point.padding = 0.3,
      # segment.color = "grey50",
      # segment.size = 0.3,
      # max.overlaps = 20,
      # ylim = c(150, NA),
      # inherit.aes = FALSE
    )
}


# Annotate the top-scoring SNP with a horizontal guide to y-axis and value label
top_snp <- META %>%
  arrange(desc(log10p)) %>%
  slice(1)

x_axis_pos <- min(META$bp_cum, na.rm = TRUE)
x_range <- diff(range(META$bp_cum, na.rm = TRUE))
label_x <- x_axis_pos - 0.01 * x_range

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


# Save the plot
# Update output paths as needed, "svg", "eps"
for (format in c("png", "pdf")) {
  # if it's png, set dpi to 350
  if (format == "png") {
    ggsave(paste0(output_prefix, ".", format), width = figwidth, height = figheight, 
           plot = manhplot, dpi = 350)
  } else {
    ggsave(paste0(output_prefix, ".", format), width = figwidth, height = figheight, 
         plot = manhplot)
  }
}

cat("Figure 1B manhattan plot saved successfully!\n")
