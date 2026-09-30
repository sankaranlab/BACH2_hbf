#!/usr/bin/env Rscript

suppressPackageStartupMessages({
	library(data.table)
	library(ggplot2)
	library(scales)
})

script_path <- sub("^--file=", "", 
                   grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
script_dir  <- if (length(script_path)) dirname(script_path) else getwd()
root <- system2("git", 
                c("-C", script_dir, "rev-parse", "--show-toplevel"), 
                stdout = TRUE)

## ---------------------------------------------------------------------------
## Inputs / outputs
## ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

input_path <- if (length(args) >= 1) {
	args[[1]]
} else {
	file.path(
	root,
	"StatGen_analysis/COJO_and_Finemap/EDFig2d_source_data_nonzero-cs.tsv"
  )
}

out_plot_path <- if (length(args) >= 2) {
	args[[2]]
} else {
	file.path(
	root,
	"StatGen_analysis/COJO_and_Finemap/EDFig2d_variant_function_barplot_noCS0.pdf"
  )
}

out_summary_path <- if (length(args) >= 3) {
	args[[3]]
} else {
	file.path(
	root,
	"StatGen_analysis/COJO_and_Finemap/EDFig2d_variant_function_summary_noCS0.tsv"
  )
}

y_mode <- if (length(args) >= 4) {
	tolower(args[[4]])
} else {
	"proportion"
}

if (!y_mode %in% c("proportion", "count")) {
	stop("4th argument (y_mode) must be either 'proportion' or 'count'")
}

## ---------------------------------------------------------------------------
## Load and process data
## ---------------------------------------------------------------------------

ignore_cs0 <- TRUE

dt <- fread(input_path)

if (!"functional_element" %in% names(dt)) {
	stop("Input file must contain a functional_element column")
}

n_removed_cs0 <- 0L
if (ignore_cs0) {
	if (!"CREDIBLE_SET" %in% names(dt)) {
		stop("ignore_cs0=TRUE requires a CREDIBLE_SET column in the input file")
	}
	dt[, CREDIBLE_SET_numeric := suppressWarnings(as.integer(CREDIBLE_SET))]
	n_before <- nrow(dt)
	dt <- dt[!is.na(CREDIBLE_SET_numeric) & CREDIBLE_SET_numeric != 0L]
	n_removed_cs0 <- n_before - nrow(dt)
	dt[, CREDIBLE_SET_numeric := NULL]

	if (nrow(dt) == 0) {
		stop("All rows were removed after filtering CREDIBLE_SET == 0")
	}
}

category_order <- c("3'UTR", "TSS", "exon", "intron", "intergenic", "promoter")

category_colors <- c(
	"3'UTR" = "#44B27A",
	"TSS" = "#F39A2B",
	"exon" = "#F2D8C7",
	"intron" = "#3E7B9B",
	"intergenic" = "#CFE8F3",
	"promoter" = "#C9965B"
)

dt[, functional_element := trimws(functional_element)]
dt[, functional_element := fifelse(
	is.na(functional_element) | functional_element == "",
	"intergenic",
	functional_element
)]

dt[, functional_element := fcase(
	tolower(functional_element) == "3'utr", "3'UTR",
	toupper(functional_element) == "TSS", "TSS",
	tolower(functional_element) == "exon", "exon",
	tolower(functional_element) == "intron", "intron",
	tolower(functional_element) == "intergenic", "intergenic",
	tolower(functional_element) == "promoter", "promoter",
	default = "intergenic"
)]

## ---------------------------------------------------------------------------
## Generate summary and create plot
## ---------------------------------------------------------------------------

summary_dt <- dt[
	, .(n_variants = .N),
	by = .(functional_element)
][
	, proportion := n_variants / sum(n_variants)
]

summary_dt[, functional_element := factor(functional_element, levels = category_order)]
setorder(summary_dt, functional_element)

summary_out <- copy(summary_dt)
summary_out[, functional_element := as.character(functional_element)]
fwrite(summary_out, out_summary_path, sep = "\t")

if (y_mode == "proportion") {
	p <- ggplot(summary_dt, aes(x = "All variants", y = proportion, fill = functional_element)) +
		geom_col(width = 0.32, color = NA) +
		scale_fill_manual(values = category_colors, drop = FALSE, name = "Variant category") +
		scale_y_continuous(
			labels = percent_format(accuracy = 1),
			expand = expansion(mult = c(0, 0.02))
		) +
		labs(x = NULL, y = "Proportion of variants") +
		theme_classic(base_size = 12) +
		theme(
			axis.text.x = element_blank(),
			axis.ticks.x = element_blank(),
			legend.title = element_text(size = 11),
			legend.text = element_text(size = 10),
			legend.position = "right"
		)
} else {
	p <- ggplot(summary_dt, aes(x = "All variants", y = n_variants, fill = functional_element)) +
		geom_col(width = 0.32, color = NA) +
		scale_fill_manual(values = category_colors, drop = FALSE, name = "Variant category") +
		scale_y_continuous(
			labels = label_number(),
			expand = expansion(mult = c(0, 0.02))
		) +
		labs(x = NULL, y = "Number of variants") +
		theme_classic(base_size = 12) +
		theme(
			axis.text.x = element_blank(),
			axis.ticks.x = element_blank(),
			legend.title = element_text(size = 11),
			legend.text = element_text(size = 10),
			legend.position = "right"
		)
}

ggsave(
	filename = out_plot_path,
	plot = p,
	width = 4.5,
	height = 6,
	# dpi = 350
)

message("Input: ", input_path)
message("Total variants counted: ", nrow(dt))
message("Y-axis mode: ", y_mode)
message("Ignore CREDIBLE_SET == 0: ", ignore_cs0, " (default)")
message("Rows removed with CREDIBLE_SET == 0: ", n_removed_cs0)
message("Summary written to: ", out_summary_path)
message("Plot written to: ", out_plot_path)
