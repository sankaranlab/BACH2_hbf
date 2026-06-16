library(tidyverse)
library(ggpubr)
library(data.table)
library(ggrepel)
library(RColorBrewer)
library(patchwork)

b36_GENE_REF <- '/path/to/hg19_ncbiRefSeq_RefSeqAll_2026UCSC.tsv.gz'  # Update this path as needed

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

# hard-code hg19 MAMA files here:
# header: 'hg19_chr', 'hg19_pos', 'hg19_id', 'rsID', 'SNP', 'hg38_chr', 'hg38_pos', 'A1', 'A2', 'FREQ', 'BETA', 'SE', 'Z', 'P', 'N_EFF', 'N_ORIG'
mama_files <- list()
mama_files[["EUR"]] <- "/path/to/BACH2_hbf/StatGen_analysis/MultiSuSiE/raw/hg19_sumstats/MAMA_EUR_SNPaligned_16865_hg19_cleaned_merged.tsv.gz"
mama_files[["AFR"]] <- "/path/to/BACH2_hbf/StatGen_analysis/MultiSuSiE/raw/hg19_sumstats/MAMA_AFR_SNPaligned_16865_hg19_cleaned_merged.tsv.gz"
mama_files[["THAI"]] <- "/path/to/BACH2_hbf/StatGen_analysis/MultiSuSiE/raw/hg19_sumstats/MAMA_THAI_SNPaligned_16865_hg19_cleaned_merged.tsv.gz"


if(! interactive()){
    # Expected CLI args:
    # 1) z_file_prefix 2) chrom 3) start 4) end 5) susie_result_prefix 6) gene_name
    # 7) [optional] comma-separated credible-set rsIDs to annotate, e.g. "rs1,rs2,rs3"
    # 8) [optional] leading SNP rsID for LD coloring, e.g., "rs123456"
    #    If omitted, determined as the variant with highest -log10p across all populations.

    args <- commandArgs(trailingOnly = TRUE)
    
    input_file_prefix <- args[1]

    chrom <- args[2]
    left <- as.numeric(args[3])
    right <- as.numeric(args[4])

    susie_result_prefix <- args[5] 
    gene_name <- args[6]
    annotate_rsids_arg <- args[7]
    leading_snp_arg <- args[8]

    if (length(args) > 9){
        plot_width <- as.numeric(args[9])
        plot_height <- as.numeric(args[10])
    } else{
        plot_width <- 5
        plot_height <- 4
    }

    # Start assembling figname
    figname <- paste0(gene_name, "_", chrom, "_", left, "_", right, "_MAMA-pvals_", plot_width, "x", plot_height, ".pdf")

    # Parse comma-separated rsID list for annotation; NULL means "default labeling behavior".
    if (!is.na(annotate_rsids_arg) && nchar(trimws(annotate_rsids_arg)) > 0) {
        annotate_rsids <- trimws(strsplit(annotate_rsids_arg, ",")[[1]])
        cat("Annotating specified rsIDs:", paste(annotate_rsids, collapse = ", "), "\n")
        # Add "_anntSNPs" to figname (before file extension) if specific rsIDs are being annotated.
        figname <- sub("(\\.pdf)$", "_anntSNPs\\1", figname)
    } else {
        annotate_rsids <- NULL
    }
    # Parse leading SNP argument; if not provided, will determine from data. If provided, add "_leadSNP" to figname.
    if (is.na(leading_snp_arg) || nchar(trimws(leading_snp_arg)) == 0) {
        leading_snp <- NULL
    } else {
        leading_snp <- trimws(leading_snp_arg)
        figname <- sub("(\\.pdf)$", paste0("_", leading_snp, "lead\\1"), figname)
        cat("Using specified leading SNP:", leading_snp, "\n")
    }

    # Load MAMA results for the region from each population and compute -log10p for plotting.
    eur_dt <- read_mama_region(mama_files[["EUR"]]) %>%
        mutate(log10p_eur = -log10(p))
    afr_dt <- read_mama_region(mama_files[["AFR"]]) %>%
        mutate(log10p_afr = -log10(p))
    thai_dt <- read_mama_region(mama_files[["THAI"]]) %>%
        mutate(log10p_thai = -log10(p))
        gwas_dt <- merge(eur_dt, afr_dt, by = c("rsid", "chr", "position"), all = TRUE) %>%
             merge(thai_dt, by = c("rsid", "chr", "position"), all = TRUE)

    # Determine leading SNP if not specified on command line.
    if (is.na(leading_snp_arg) || nchar(trimws(leading_snp_arg)) == 0) {
        leading_snp <- find_leading_snp(gwas_dt)
        print(sprintf("No leading SNP specified; using %s (highest -log10p = %.2f)", 
                      leading_snp, 
                      max(gwas_dt$log10p_eur, gwas_dt$log10p_afr, gwas_dt$log10p_thai, na.rm = TRUE)))
    } else {
        leading_snp <- trimws(leading_snp_arg)
        print(sprintf("Using specified leading SNP: %s", leading_snp))
    }

    # Extract LD scores for each population.
    # need z files to get variant names and positions in the same order as the LD matrix
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
        annotate_positions <- sort(unique(gwas_dt$position[gwas_dt$rsid %in% annotate_rsids]))
        if (length(annotate_positions) == 0) {
            warning("None of the requested annotate_rsids were found in GWAS data for this region.")
        }
    } else {
        annotate_positions <- NULL
    }
     
    # SuSiE outputs 
    credible_set_file <- paste0(susie_result_prefix, "_credible_sets.tsv") # header: component       num_variants    purity  coverage        passed_filter   variant_indices variants
    
    # Use a shared set of x-axis breaks so vertical gridlines align across all tracks.
    x_breaks <- pretty(c(left, right), n = 5)

    p_gwas <- plot_gwas_tracks(gwas_dt, chrom, left, right, 
                               credible_set_file,
                               x_breaks = x_breaks, annotate_rsids = annotate_rsids,
                               annotate_positions = annotate_positions, ld_data_list = ld_data_list,
                               leading_snp = leading_snp, mode = "mama") #+ ggtitle(gene_name)

    
    # Build gene structures in the same coordinate window as association tracks.
    genes <- getGenes(chrom, left, right)
    ## adjust dataframe for plotting
    genes = genes[, .(start=pmax(txStart, left), 
                      end=pmin(txEnd, right), 
                      cdsStart, cdsEnd, name=name2, 
                      exonStarts, exonEnds)]
    print(dim(genes))

    if(nrow(genes) == 0){
        # report error and quit if no gene is found in the region
        stop("No gene found in the specified region.")
    }

    # Plot transcript/exon track with non-overlapping label layers.
    print(unique(genes$name))
    print('Plotting genes...')
    pgList <- plotFancyGenes(genes)
    pg <- pgList$plot
    layers <- pgList$nLayer
    ## add aesthetics
    pg <- pg + scale_x_continuous(limits = c(left, right), breaks = x_breaks,
                                    labels = function(x) x / 1e6,
                                    expand = expansion(mult = c(0.01, 0.01))) + 
        theme(panel.grid.major.x=element_line(linewidth=0.3,color='#cccccc'),
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
    ## remove x axis for gene track
    pg <- pg  + theme(axis.title.y=element_blank(),
                      axis.text.y=element_blank(), 
                      # axis.line.x=element_blank(),
                      axis.title.x=element_blank(), 
                      axis.text.x=element_blank(), 
                      axis.line.x=element_blank(),
                      axis.line.y=element_line(color='black',linewidth=0.5)
                    )
    ## remove y axis texts and labels and add x axis for gene plot
    p_gwas <- p_gwas + xlab(paste('(hg19) Position on Chromosome', chrom, '(Mb)')) + 
            theme(axis.ticks=element_blank(), 
                #   axis.line.y=element_line(color='black',linewidth=0.5),
                  axis.title.x=element_text(size=10,hjust=0.5), 
                  axis.text.x=element_text(size=8,color='black'),
                  axis.line=element_line(color='black',linewidth=0.5), 
                  plot.margin=margin(b=5,unit='pt'),
                #   plot.margin=margin(t=-1,r=25,b=0,l=35,unit='pt')
                  legend.title=element_text(size=8),
                  legend.text=element_text(size=7)
                  )

        if (!is.null(annotate_positions) && length(annotate_positions) > 0) {
                pg <- pg +
                        geom_vline(xintercept = annotate_positions, linetype = "dashed", color = "#2b8cbe", linewidth = 0.35, alpha = 0.7)
        }

    # Combine GWAS and gene tracks.
    # patchwork aligns inner panel widths correctly regardless of y-axis label width.
    p_combined <- wrap_plots(c(list(pg, p_gwas)), ncol = 1,
                             heights = c(0.15, 3))

    # Save to different file types
    for (ext in c("pdf", "png")) { #, "svg", "tiff"
        out_file <- sub("\\.pdf$", paste0(".", ext), figname)
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