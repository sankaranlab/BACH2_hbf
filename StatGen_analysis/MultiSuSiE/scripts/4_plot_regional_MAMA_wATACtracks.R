library(tidyverse)
library(ggpubr)
library(data.table)
library(ggrepel)
library(RColorBrewer)
library(rtracklayer)
library(patchwork)

b36_GENE_REF <- '/path/to/hg19_ncbiRefSeq_RefSeqAll_2026UCSC.tsv.gz'  # Update this path as needed

ATAC_DIR_WENG2024   <- '/path/to/bone_marrow_atac_peak_bigwig/hg19_liftOver/'  # Update this path as needed
ATAC_DIR_CORCES2016 <- '/path/to/Corces_2016_hg19bw/'  # Update this path as needed

# Load shared step4 utility functions 
# source("scripts/step4_utility.R")
get_script_dir <- function() {
    file_arg <- sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))])
    if (length(file_arg) > 0) {
        return(dirname(normalizePath(file_arg[1])))
    }
    return(getwd())
}
source(file.path(get_script_dir(), "step4_utility.R"))

# Hard-code MAMA files here:
# header: 'hg19_chr', 'hg19_pos', 'hg19_id', 'rsID', 'SNP', 'hg38_chr', 'hg38_pos', 'A1', 'A2', 'FREQ', 'BETA', 'SE', 'Z', 'P', 'N_EFF', 'N_ORIG'
mama_files <- list()
mama_files[["EUR"]] <- "/path/to/BACH2_hbf/StatGen_analysis/MultiSuSiE/raw/hg19_sumstats/MAMA_EUR_SNPaligned_16865_hg19_cleaned_merged.tsv.gz"
mama_files[["AFR"]] <- "/path/to/BACH2_hbf/StatGen_analysis/MultiSuSiE/raw/hg19_sumstats/MAMA_AFR_SNPaligned_16865_hg19_cleaned_merged.tsv.gz"
mama_files[["THAI"]] <- "/path/to/BACH2_hbf/StatGen_analysis/MultiSuSiE/raw/hg19_sumstats/MAMA_THAI_SNPaligned_16865_hg19_cleaned_merged.tsv.gz"


if(! interactive()){
    # Expected CLI args:
    # 1) input_file_prefix 2) chrom 3) start 4) end 5) susie_result_prefix 6) gene_name
    # 7) [optional] comma-separated cell types to plot, e.g. "CD14,CD34,HSC"
    #    If omitted or empty, all BigWig/bedgraph files in the chosen ATAC directory are plotted.
    # 8) [optional] comma-separated credible-set rsIDs to annotate, e.g. "rs1,rs2,rs3"
    #    If omitted or empty, default labeling behavior is used.
    # 9) [optional] leading SNP rsID for LD coloring, e.g., "rs123456"
    #    If omitted, determined as the variant with highest -log10p across EUR/AFR/THAI.
    # 10) [optional] ATAC data source: "Corces2016" (default) or "Weng2024"
    input_file_prefix <- commandArgs(trailingOnly = TRUE)[1]

    chrom <- commandArgs(trailingOnly = TRUE)[2]
    left <- as.numeric(commandArgs(trailingOnly = TRUE)[3])
    right <- as.numeric(commandArgs(trailingOnly = TRUE)[4])

    susie_result_prefix <- commandArgs(trailingOnly = TRUE)[5] 
    gene_name <- commandArgs(trailingOnly = TRUE)[6]

    cell_type_filter_arg <- commandArgs(trailingOnly = TRUE)[7]

    annotate_rsids_arg <- commandArgs(trailingOnly = TRUE)[8]
    leading_snp_arg    <- commandArgs(trailingOnly = TRUE)[9]
    atac_source_arg    <- commandArgs(trailingOnly = TRUE)[10]  # "Corces2016" (default) or "Weng2024"

    # Set ATAC data directory and figure label based on source selection.
    print(sprintf("Selected ATAC source: %s", ifelse(is.na(atac_source_arg), "Corces2016 (default)", atac_source_arg)))
    atac_source <- if (is.na(atac_source_arg)) "Corces2016" else trimws(atac_source_arg)
    if (atac_source == "Weng2024") {
        bedgraph_path <- ATAC_DIR_WENG2024
        atac_label    <- "Weng2024tracks"
    } else if (atac_source == "Corces2016") {
        bedgraph_path <- ATAC_DIR_CORCES2016
        atac_label    <- "Corces2016tracks"
    } else {
        stop("Invalid ATAC source specified. Use 'Corces2016' or 'Weng2024'.")
    }
    print(atac_label)

    figname <- paste0(gene_name, "_MAMA-pval_", chrom, "_", left, "_", right, "_combined_", atac_label, ".png")

    # Parse comma-separated cell type list; NULL means "plot all".
    if (!is.na(cell_type_filter_arg) && 
        nchar(trimws(cell_type_filter_arg)) > 3) {
        cell_type_filter <- trimws(strsplit(cell_type_filter_arg, ",")[[1]])
        # count number of cell types for figure naming
        num_cell_types <- length(cell_type_filter)
        figname <- paste0(gene_name, "_MAMA-pval_", chrom, "_", left, "_", right, "_combined_", num_cell_types, "_", atac_label, ".png")
    } else {
        cell_type_filter <- NULL
    }
  
    # Parse comma-separated rsID list for annotation; NULL means "default labeling behavior".
    if (!is.na(annotate_rsids_arg) && nchar(trimws(annotate_rsids_arg)) > 3) {
        annotate_rsids <- trimws(strsplit(annotate_rsids_arg, ",")[[1]])
        # Add "_anntSNPs" to figname (before file extension) if specific rsIDs are being annotated.
        figname <- sub("(\\.png)$", "_anntSNPs\\1", figname)
    } else {
        annotate_rsids <- NULL
    }

    eur_dt <- read_mama_region(mama_files[["EUR"]]) %>%
        mutate(log10p_eur = -log10(p))
    afr_dt <- read_mama_region(mama_files[["AFR"]]) %>%
        mutate(log10p_afr = -log10(p))
    thai_dt <- read_mama_region(mama_files[["THAI"]]) %>%
        mutate(log10p_thai = -log10(p))
    gwas_dt <- merge(eur_dt, afr_dt, by = c("rsid", "chr", "position"), all = TRUE) %>%
               merge(thai_dt, by = c("rsid", "chr", "position"), all = TRUE)

    # Determine leading SNP if not specified on command line.
    if (is.na(leading_snp_arg) || nchar(trimws(leading_snp_arg)) <= 3) {
        leading_snp <- find_leading_snp(gwas_dt)
        print(sprintf("No leading SNP specified; using %s (highest -log10p = %.2f)", 
                      leading_snp, 
                      max(gwas_dt$log10p_eur, gwas_dt$log10p_afr, gwas_dt$log10p_thai, na.rm = TRUE)))
    } else {
        leading_snp <- trimws(leading_snp_arg)
        print(sprintf("Using specified leading SNP: %s", leading_snp))
    }

    # Extract LD scores for each population.
    # Use z files to keep rsID order consistent with LD matrices.
    ld_data_list <- list()
    for (pop in c("EUR", "AFR", "THAI")) {
        z_file <- paste0(input_file_prefix, ".", pop, ".z")
        ld_file <- paste0(input_file_prefix, ".", pop, ".ld")
        if (file.exists(z_file) && file.exists(ld_file)) {
            pop_dt <- fread(z_file) %>%
                select(rsid, position, p) %>%
                mutate(log10p = -log10(p))
            ld_data_list[[pop]] <- extract_ld_stats(ld_file, pop_dt, leading_snp)
        }
    }

    # Resolve designated rsIDs to genomic positions once and reuse for all panels.
    if (!is.null(annotate_rsids)) {
        cat("Annotating specified rsIDs: ", paste(annotate_rsids, collapse = ", "), "\n")
        annotate_positions <- sort(unique(gwas_dt$position[gwas_dt$rsid %in% annotate_rsids]))
        if (length(annotate_positions) == 0) {
            warning("None of the requested annotate_rsids were found in GWAS data for this region.")
        }
    } else {
        cat("No specific rsIDs provided for annotation; using default labeling behavior.\n")
        annotate_positions <- NULL
    }
     
    # SuSiE outputs used to annotate credible-set variants.
    credible_set_file <- paste0(susie_result_prefix, "_credible_sets.tsv") # header: component       num_variants    purity  coverage        passed_filter   variant_indices variants
    
    # Use a shared set of x-axis breaks so vertical gridlines align across all tracks.
    x_breaks <- pretty(c(left, right), n = 5)

    p_gwas <- plot_gwas_tracks(gwas_dt, chrom, left, right, credible_set_file,
                               x_breaks = x_breaks, annotate_rsids = annotate_rsids,
                               annotate_positions = annotate_positions,
                               ld_data_list = ld_data_list,
                               leading_snp = leading_snp,
                               mode = "mama") #+ ggtitle(gene_name)

    
    # Build gene structures in the same coordinate window as association tracks.
    genes <- getGenes(chrom, left, right)
    ## adjust dataframe for plotting
    genes = genes[, .(start=pmax(txStart, left), end=pmin(txEnd, right), cdsStart, cdsEnd, name=name2, exonStarts, exonEnds)]
    print(dim(genes))
    if(nrow(genes) == 0){
        # report error and quit if no gene is found in the region
        stop("No gene found in the specified region.")
    }

    # Plot transcript/exon track with non-overlapping label layers.
    print(unique(genes$name))
    print('Plotting genes...')
    pgList <- plotFancyGenes(genes,
                             label_force = 25,
                             label_size = 4,
                             label_family = 'Helvetica',
                             label_ylim = c(NA, -6))
    pg <- pgList$plot
    layers <- pgList$nLayer
    ## add aesthetics
    pg <- pg + scale_x_continuous(limits = c(left, right), breaks = x_breaks, expand = expansion(mult = c(0.01, 0.01))) + 
        theme(panel.grid.major.x=element_line(linewidth=0.3, color='#cccccc'),
              panel.grid.minor=element_blank(),
              panel.grid.major.y=element_blank(),
              axis.line.y=element_line(color='black', linewidth=0.5),
              plot.margin=margin(t=0,r=25,b=5,l=15,unit='pt')
              )
    if (layers > 10){
        n=3
    }else if (layers>5){
        n = 2.5
    }else{
        n=0.5*layers
    }
    ## remove x axis for score plot
    pg <- pg  + theme(axis.title.y=element_blank(),
                  axis.text.y=element_blank(), # axis.line.x=element_blank(),
                  axis.ticks=element_blank(), 
                  axis.title.x=element_blank(), 
                  axis.text.x=element_blank(), 
                  axis.line.x=element_blank(),
                  axis.line.y=element_line(color='black', linewidth=0.5)
                 )
    ## remove y axis texts and labels and add x axis for gwas plot
    p_gwas <- p_gwas + xlab(paste('(hg19) Position on Chromosome',chrom) ) + 
            theme(
                #   axis.line.y=element_line(color='black', linewidth=0.5),
                  axis.title.x=element_text(size=12,hjust=0.5), 
                  axis.text.x=element_text(size=10,color='black'),
                  axis.line=element_line(color='black', linewidth=0.5), 
                  plot.margin=margin(b=5,unit='pt'),
                  legend.title=element_text(size=9), legend.text=element_text(size=8)
                #   plot.margin=margin(t=-1,r=25,b=0,l=35,unit='pt')
                  )

        if (!is.null(annotate_positions) && length(annotate_positions) > 0) {
                pg <- pg +
                        geom_vline(xintercept = annotate_positions, linetype = "dashed", color = "#2b8cbe", linewidth = 0.35, alpha = 0.7)
        }

    # Plot each ATAC track (BigWig or bedgraph) as an additional vertical track.
    atac_files <- list.files(bedgraph_path,
                             pattern = "\\.(bw|bedgraph)$",
                             full.names = TRUE, ignore.case = TRUE)
    # Optionally restrict to a subset of cell types (matched by display label).
    if (!is.null(cell_type_filter)) {
        atac_labels <- get_cell_type(atac_files)
        keep_idx <- which(atac_labels %in% cell_type_filter)
        # Order selected tracks to follow the user-specified cell_type_filter order.
        ordered_keep_idx <- keep_idx[order(match(atac_labels[keep_idx], cell_type_filter), keep_idx)]
        filtered_atac_files <- atac_files[ordered_keep_idx]
        if (length(filtered_atac_files) == 0) {
            print(get_cell_type(atac_files))
            stop("No ATAC files matched the requested cell types: ", paste(cell_type_filter, collapse = ", "))
        }
        atac_files <- filtered_atac_files
    }
    # count the number of files
    num_beds <- length(atac_files)
    # iterate through the file list, save plots
    bg_tracks <- list()
    for (atac_file in atac_files) {
        p_bg <- plot_bigwig_track(atac_file, chrom, left, right, x_breaks = x_breaks,
                      annotate_positions = annotate_positions,
                      bar_fill = "gray10")
        # aesthetics
        p_bg <- p_bg + theme(
            axis.title.x = element_blank(),
            axis.text.y = element_blank(),
            panel.grid.major.y = element_blank(),
            # plot.margin = margin(t=0, r=25, b=5, l=15, unit='pt')
        )
        bg_tracks <- c(bg_tracks, list(p_bg))
    }
    
    # Stack GWAS panel on top of all signal tracks.
    # Extract legend and place it above the gene track.
    legend <- get_legend(p_gwas)
    p_gwas_no_legend <- p_gwas + theme(legend.position = "none")
    # patchwork aligns inner panel widths correctly regardless of y-axis label width.
    p_combined <- wrap_plots(c(list(legend, pg, p_gwas_no_legend), bg_tracks), ncol = 1,
                             heights = c(0.5, 0.5, 4, rep(0.5, num_beds)))
    
    plot_width <- 5
    plot_height <- 4 + num_beds * 0.5  # Base height plus additional
    # Save to different file types
    for (ext in c("pdf", "png")) { #, "svg", "tiff"
        out_file <- sub("\\.png$", paste0(".", ext), figname)
        cat(sprintf("Saving plot to %s...\n", out_file))
        if (ext == "png") {
            ggsave(out_file, p_combined, 
                   width = plot_width, height = plot_height, 
                   units = "in", dpi = 350) # , limitsize = FALSE
        } else {
            ggsave(out_file, p_combined, 
                width = plot_width, height = plot_height, 
                units = "in") # , limitsize = FALSE, dpi = 350
        }
    }
}