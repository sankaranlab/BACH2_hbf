#!/usr/bin/env Rscript

suppressPackageStartupMessages({
	library(data.table)
	library(ggplot2)
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
	"StatGen_analysis/COJO_and_Finemap/EDFig2c_source_data_5e-8.tsv"
  )
}

out_plot_path <- if (length(args) >= 2) {
	args[[2]]
} else {
	file.path(
	root,
	"StatGen_analysis/COJO_and_Finemap/EDFig2c_varcount_by_gene_barplot.pdf"
  )
}

out_summary_path <- if (length(args) >= 3) {
	args[[3]]
} else {
	file.path(
	root,
	"StatGen_analysis/COJO_and_Finemap/EDFig2c_varcount_by_gene_summary.tsv"
  )
}

dt <- fread(input_path)

required_cols <- c("window_name", "fema_cojo_ncausal", "finemap_credset")
missing_cols <- setdiff(required_cols, names(dt))
if (length(missing_cols) > 0) {
	stop("Input is missing required columns: ", paste(missing_cols, collapse = ", "))
}

## ---------------------------------------------------------------------------
## Prepare data for plotting
## ---------------------------------------------------------------------------

plot_dt <- dt[, .(
	window_name,
	fema_cojo_ncausal = as.integer(fema_cojo_ncausal),
	finemap_credset = as.integer(finemap_credset)
)]

## Use the first symbol when multiple genes are joined (e.g., CTC1;KRBA2 -> CTC1)
plot_dt[, gene := sub(";.*$", "", window_name)]

long_dt <- melt(
	plot_dt,
	id.vars = c("window_name", "gene"),
	measure.vars = c("fema_cojo_ncausal", "finemap_credset"),
	variable.name = "metric",
	value.name = "n_variants"
)

## ---------------------------------------------------------------------------
## Format and create visualization
## ---------------------------------------------------------------------------

metric_labels <- c(
	"fema_cojo_ncausal" = "COJO",
	"finemap_credset" = "95% Credible Set"
)

long_dt[, metric := factor(metric, levels = names(metric_labels))]
setorder(long_dt, metric, -n_variants, gene)

## Reorder genes independently within each panel
long_dt[, gene_in_panel := paste(metric, gene, sep = "___")]
long_dt[, gene_in_panel := factor(gene_in_panel, levels = unique(gene_in_panel))]

summary_dt <- dcast(
	long_dt,
	gene ~ metric,
	value.var = "n_variants"
)
setnames(summary_dt, old = names(metric_labels), new = unname(metric_labels))
setorder(summary_dt, -`COJO`, -`95% Credible Set`, gene)
fwrite(summary_dt, out_summary_path, sep = "\t")

p <- ggplot(long_dt, aes(x = n_variants, y = gene_in_panel)) +
	geom_col(width = 0.85, fill = "#666666", color = "#666666") +
	facet_wrap(~metric, scales = "free", nrow = 1, labeller = as_labeller(metric_labels)) +
	scale_y_discrete(labels = function(x) sub("^.*___", "", x)) +
	labs(x = "Number of variants", y = NULL) +
	theme_classic(base_size = 10) +
	theme(
		strip.background = element_blank(),
		strip.text = element_text(size = 8),
		panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
		panel.spacing.x = unit(0.9, "lines"),
		axis.text.y = element_text(size = 10),
		axis.text.x = element_text(size = 10),
		axis.title.x = element_text(size = 10)
	)

ggsave(
	filename = out_plot_path,
	plot = p,
	width = 5,
	height = 5,
	dpi = 350
)

message("Input: ", input_path)
message("Genes plotted: ", uniqueN(plot_dt$gene))
message("Summary written to: ", out_summary_path)
message("Plot written to: ", out_plot_path)
