# Merge EUR and AFR COJO results, get Z score fold change,
# filter out variants with inflated ZJ compared to Z, and annotate nearest genes.
library(data.table)
library(tidyverse)

# ref header: bin	name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	score	name2	cdsStartStat	cdsEndStat	exonFrames
# Name is the RefSeq ID, name2 is the gene symbol
ref_genes_path <- '/lab-share/Hem-Sankaran-e2/Public/ref_genomes/human/hg38_ncbiRefSeq_AllGenes_2026UCSC.tsv.gz'

normalize_chr <- function(x) {
    chr <- as.character(x)
    chr <- sub('^chr', '', chr, ignore.case = TRUE)
    toupper(chr)
}

get_z_and_foldchange <- function(cojo_dt) {
    # Read in the dt with at least the following columns: b, se, bJ, bJ_se
    # z := b/se, zJ := bJ/bJ_se, zJfc := abs(zJ/z)
    # return z, zJ, and zJfc

    cojo_dt <- as.data.table(copy(cojo_dt))
    required_cols <- c('b', 'se', 'bJ', 'bJ_se')
    missing_cols <- setdiff(required_cols, names(cojo_dt))
    if (length(missing_cols) > 0) {
        stop('COJO table is missing required columns: ', paste(missing_cols, collapse = ', '))
    }

    cojo_dt[, z := fifelse(is.na(se) | se == 0, NA_real_, b / se)]
    cojo_dt[, zJ := fifelse(is.na(bJ_se) | bJ_se == 0, NA_real_, bJ / bJ_se)]
    cojo_dt[, zJfc := fifelse(is.na(z) | z == 0 | is.na(zJ), NA_real_, abs(zJ / z))]

    cojo_dt
}

annotate_nearest_genes <- function(merged_df, ref_genes_path) {
    # Annotate nearest genes for the merged df using the reference genes.
    # Return the annotated df.

    merged_dt <- as.data.table(copy(merged_df))
    if (!'Chr' %in% names(merged_dt)) {
        stop('Input table for gene annotation must contain a Chr column')
    }

    pos_col <- if ('bp' %in% names(merged_dt)) {
        'bp'
    } else if ('Position' %in% names(merged_dt)) {
        'Position'
    } else {
        stop('Input table for gene annotation must contain either bp or Position')
    }

    if (nrow(merged_dt) == 0) {
        return(character())
    }

    merged_dt[, row_id := .I]
    merged_dt[, chr_norm := normalize_chr(Chr)]
    merged_dt[, pos_bp := suppressWarnings(as.integer(get(pos_col)))]

    ref_dt <- fread(
        ref_genes_path,
        select = c('name', 'name2', 'chrom', 'txStart', 'txEnd')
    )
    ref_dt <- ref_dt[
        grepl('^NM', name) & !is.na(chrom) & !is.na(txStart) & !is.na(txEnd),
        .(
            refseq = name,
            gene = fifelse(is.na(name2) | name2 == '', name, name2),
            ref_chr = normalize_chr(chrom),
            tx_start = as.integer(txStart) + 1L,
            tx_end = as.integer(txEnd)
        )
    ]

    if (nrow(ref_dt) == 0) {
        stop('No NM RefSeq entries found in reference: ', ref_genes_path)
    }

    query_dt <- merged_dt[!is.na(chr_norm) & !is.na(pos_bp)]
    nearest_gene <- rep(NA_character_, nrow(merged_dt))
    if (nrow(query_dt) == 0) {
        return(nearest_gene)
    }

    inside_hits <- ref_dt[
        query_dt,
        on = .(ref_chr = chr_norm, tx_start <= pos_bp, tx_end >= pos_bp),
        nomatch = 0L,
        allow.cartesian = TRUE
    ]

    inside_best <- inside_hits[
        order(row_id, gene, refseq),
        .SD[1],
        by = row_id
    ][
        , .(row_id, nearest_gene = gene)
    ]

    unmatched <- query_dt[!inside_best, on = 'row_id']

    ref_end <- copy(ref_dt)
    setorder(ref_end, ref_chr, tx_end, gene, refseq)
    setkey(ref_end, ref_chr, tx_end)

    ref_start <- copy(ref_dt)
    setorder(ref_start, ref_chr, tx_start, gene, refseq)
    setkey(ref_start, ref_chr, tx_start)

    left_hit <- ref_end[
        unmatched,
        on = .(ref_chr = chr_norm, tx_end = pos_bp),
        roll = Inf
    ]

    right_hit <- ref_start[
        unmatched,
        on = .(ref_chr = chr_norm, tx_start = pos_bp),
        roll = -Inf
    ]

    dist_left <- ifelse(is.na(left_hit$tx_end), Inf, unmatched$pos_bp - left_hit$tx_end)
    dist_right <- ifelse(is.na(right_hit$tx_start), Inf, right_hit$tx_start - unmatched$pos_bp)
    left_key <- paste0(left_hit$gene, '|', left_hit$refseq)
    right_key <- paste0(right_hit$gene, '|', right_hit$refseq)

    choose_left <- dist_left < dist_right
    choose_right <- dist_right < dist_left
    tie <- is.finite(dist_left) & is.finite(dist_right) & (dist_left == dist_right)
    if (any(tie)) {
        choose_left[tie] <- left_key[tie] <= right_key[tie]
        choose_right[tie] <- !choose_left[tie]
    }

    unmatched_best <- unmatched[, .(row_id)]
    unmatched_best[, nearest_gene := fifelse(
        choose_left,
        left_hit$gene,
        fifelse(choose_right, right_hit$gene, NA_character_)
    )]

    annot_dt <- rbindlist(list(inside_best, unmatched_best), use.names = TRUE, fill = TRUE)
    nearest_gene[annot_dt$row_id] <- annot_dt$nearest_gene
    return(nearest_gene)
}

if(!interactive()){
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) < 4) {
        stop('Usage: Rscript merge_eur-afr_cojos.R <eur_cojo.tsv> <afr_cojo.tsv> <output_prefix> <cutoff> [ref_genes.tsv.gz]')
    }
    print(args)
    eur_cojo_path <- args[[1]]
    afr_cojo_path <- args[[2]]
    output_prefix <- args[[3]]
    cutoff <- as.numeric(args[[4]])
    if (length(args) >= 5) {
        ref_genes_path <- args[[5]]
    }
    
    cat("Using fold-change cutoff of", cutoff, "\n")

    # Header: Chr     SNP     bp      refA    freq    b       se      p       n       freq_geno       bJ      bJ_se       pJ      LD_r
    eur_dt <- fread(eur_cojo_path) %>% get_z_and_foldchange()
    afr_dt <- fread(afr_cojo_path) %>% get_z_and_foldchange()

    merged_dt <- merge(eur_dt, afr_dt, all = TRUE,
                    by = c('Chr', 'bp', 'SNP', 'refA', 'freq', 'b', 'se', 'p'), 
                    suffixes = c('_eur', '_afr'))
    print(summary(merged_dt[, .(zJfc_eur, zJfc_afr)]))

    # filter
    cat("Filtering merged results with fold-change cutoff of", cutoff, "...\n")
    cat("Total rows before filtering:", nrow(merged_dt), "\n")
    cat("\tEUR:", sum(!is.na(merged_dt$n_eur)), "non-NA rows\n")
    cat("\tAFR:", sum(!is.na(merged_dt$n_afr)), "non-NA rows\n")
    merged_df <- merged_dt %>% 
        filter(zJfc_eur <= cutoff | zJfc_afr <= cutoff) 
    cat("Total rows after filtering:", nrow(merged_df), "\n")
    cat("\tEUR:", sum(!is.na(merged_df$n_eur)), "non-NA rows\n")
    cat("\tAFR:", sum(!is.na(merged_df$n_afr)), "non-NA rows\n")
    # get nearest genes
    merged_df$Nearest_coding_gene <- annotate_nearest_genes(merged_df, ref_genes_path)
    # reorder columns
    merged_df <- merged_df %>%
        select(Chr, Nearest_coding_gene, Position=bp, rsID=SNP, 
                refA, freq, beta=b, se, z, p,
                n_eur, freq_geno_eur, bJ_eur, bJ_se_eur, zJ_eur, pJ_eur, LD_r_eur,
                n_afr, freq_geno_afr, bJ_afr, bJ_se_afr, zJ_afr, pJ_afr, LD_r_afr
                )
    # write out
    out_path <- paste0(output_prefix,"_FCcutoff", cutoff, "_annotated.tsv")
    cat("Writing merged and annotated results to", out_path, "\n")
    fwrite(merged_df, out_path, sep = "\t")
}
