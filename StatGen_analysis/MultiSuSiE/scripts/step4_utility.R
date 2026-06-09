library(data.table)
library(tidyverse)
library(ggrepel)
# Shared utility functions for step4 plotting scripts.

# Normalize chromosome labels so values like "1" and "chr1" are treated identically.
normalize_chr <- function(x) {
	x_chr <- as.character(x)
	x_chr <- trimws(x_chr)
	x_chr <- gsub('^chr', '', x_chr, ignore.case = TRUE)
	return(x_chr)
}

# Read MAMA results for the specified region from a given file path, returning a data.table with rsid, chr, position, and p-value.
read_mama_region <- function(file_path) {
    dt <- fread(file_path)
    cat("MAMA file columns:", paste(colnames(dt), collapse = ", "), ".\n")
    # Handle potential UTF-8 BOM in the first header cell.
    colnames(dt) <- sub("^\\uFEFF", "", colnames(dt))

    required_cols <- c("rsID", "hg19_chr", "hg19_pos", "P")
    missing_cols <- setdiff(required_cols, colnames(dt))
    if (length(missing_cols) > 0) {
        stop(sprintf(
            "MAMA file %s is missing required columns: %s",
            file_path,
            paste(missing_cols, collapse = ", ")
        ))
    }

    dt %>%
        select(rsID, hg19_chr, hg19_pos, P) %>%
        filter(normalize_chr(hg19_chr) == normalize_chr(chrom),
                hg19_pos >= left,
                hg19_pos <= right) %>%
        transmute(rsid = rsID,
                    chr = hg19_chr,
                    position = hg19_pos,
                    p = as.numeric(P))
}

# read gene information
getGenes <- function(chr, left, right,
					 b37_gene_ref = '/lab-share/Hem-Sankaran-e2/Public/ref_genomes/human/hg19_ncbiRefSeq_RefSeqAll_2026UCSC.tsv.gz') {
	# ref header: #bin	name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	score	name2	cdsStartStat	cdsEndStat	exonFrames
	genes <- data.table(read.table(b37_gene_ref, sep = '\t', header = TRUE))
	colnames(genes) <- c('bin', 'name', 'chrom', 'strand',
						 'txStart', 'txEnd', 'cdsStart', 'cdsEnd',
						 'exonCount', 'exonStarts', 'exonEnds', 'score',
						 'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames')
	print(dim(genes))
	chr_norm <- normalize_chr(chr)
	return(genes[normalize_chr(chrom) == chr_norm & txStart <= right & txEnd >= left])
}


# plot gene tracks
plotFancyGenes <- function(genes,
						   label_force = 60,
						   label_size = 3.5,
						   label_family = NULL,
						   label_ylim = NULL) {
  dups <- genes[, .N, by = name]
  print(dups)
  print(dups[order(-N)])
  # if any gene in the region has >3 transcript, only show the longest one
  if (max(dups$N) > 3) {
	# make each gene name unique
	genes <- genes[, .(start = start[which.max(end - start)],
					   end = end[which.max(end - start)],
					   cdStart = cdsStart[which.max(end - start)],
					   cdEnd = cdsEnd[which.max(end - start)],
					   exonStarts = exonStarts[which.max(end - start)],
					   exonEnds = exonEnds[which.max(end - start)]), by = name]
  }

  # sort the dataframe by starting position
  genes <- genes[order(start)]
  # initialize parsing
  ends <- c(0)
  Exons <- data.frame(layer = integer(), Estart = numeric(), Eend = numeric())
  Tracks <- data.frame(layer = integer(), start = numeric(), end = numeric())
  Names <- data.frame(layer = integer(), pos = numeric(), name = character())
  # check to see if two transcripts are too close and parse the exons
  for (i in 1:nrow(genes)) {
	start <- genes$start[i]
	end <- genes$end[i]
	# determine the layer to lay the gene track
	layer <- 0
	for (j in 1:length(ends)) {
	  if (ends[j] + 1000 <= start) {
		layer <- j
		ends[j] <- end
		break
	  }
	}
	if (layer == 0) {
	  ends <- c(ends, end)
	  layer <- length(ends)
	}
	track <- data.frame(layer = layer, start = start, end = end)
	Name <- data.frame(layer = layer, pos = 0.5 * (start + end), name = genes$name[i])
	Tracks <- rbind(Tracks, track)
	Names <- rbind(Names, Name)
	# parse the exon start and end positions
	exonS <- as.integer(unlist(strsplit(as.character(genes$exonStarts[i]), '[,]')))
	exonE <- as.integer(unlist(strsplit(as.character(genes$exonEnds[i]), '[,]')))

	exons <- data.frame(layer = rep(layer, length(exonS)), Estart = exonS, Eend = exonE)

	Exons <- rbind(Exons, exons)
  }
  # start plotting
  p <- ggplot() + geom_rect(data = Tracks, aes(xmin = start, xmax = end, ymin = -10 * layer + 5.5, ymax = -10 * layer + 6.5), fill = '#666666')
  p <- p + geom_rect(data = Exons, aes(xmin = Estart, xmax = Eend, ymin = -10 * layer + 4, ymax = -10 * layer + 8), fill = '#323232')

  # One label per gene: keep the first (longest) transcript row per name.
  Names_unique <- Names[!duplicated(Names$name), ]
  text_args <- list(
	data = Names_unique,
	mapping = aes(y = -10 * layer + 4, x = pos, label = name),
	force = label_force,
	size = label_size,
	fontface = 'bold.italic',
	color = 'black',
	segment.size = 0.3,
	direction = 'y'
  )
  if (!is.null(label_family)) {
	text_args$family <- label_family
  }
  if (!is.null(label_ylim)) {
	text_args$ylim <- label_ylim
  }
  p <- p + do.call(geom_text_repel, text_args)

  # final touch
  p <- p + theme_minimal() +
		   theme(axis.title.y = element_blank(),
						  axis.text.y = element_blank(),
						  axis.ticks = element_blank(),
						  panel.grid.major.x = element_line(linewidth = 0.3, color = '#cccccc'),
						  panel.grid.minor = element_blank(),
						  panel.grid.major.y = element_blank(),
						  axis.line.y = element_line(color = 'black', linewidth = 0.5)
	)

  res_list <- list(nLayer = max(Tracks$layer), plot = p)
  return(res_list)
}


# Extract LD scores for a leading SNP from an LD matrix.
# ld_file: square matrix with rsids as colnames and rownames (in same order as z_file)
# z_dt: data.table with rsid and position columns (same order as LD matrix)
# leading_rsid: the SNP to extract LD scores for
extract_ld_stats <- function(ld_file, z_dt, leading_rsid, make_symmetric = TRUE) {
	if (!file.exists(ld_file)) {
		return(NULL)
	}

	tryCatch({
		ld_raw <- fread(ld_file, header = FALSE)
		# Make it symmetric
		if (isTRUE(make_symmetric)) {
			ld_raw[lower.tri(ld_raw)] <- t(ld_raw)[lower.tri(ld_raw)]
		}

		ld_rsids <- z_dt$rsid
		ld_matrix <- as.matrix(ld_raw)

		rownames(ld_matrix) <- ld_rsids
		colnames(ld_matrix) <- ld_rsids

		# Find index of leading SNP
		lead_idx <- which(ld_rsids == leading_rsid)
		if (length(lead_idx) == 0) {
			warning(sprintf("Leading SNP %s not found in LD matrix", leading_rsid))
			return(NULL)
		}

		ld_stats <- ld_matrix[lead_idx, ]
		data.table(rsid = ld_rsids, ld_stat = as.numeric(ld_stats))
	}, error = function(e) {
		warning(sprintf("Error reading LD file %s: %s", ld_file, e$message))
		return(NULL)
	})
}


# Determine the leading SNP: highest -log10p across all populations
find_leading_snp <- function(gwas_dt, log10p_prefix = "log10p_") {
	gwas_dt_long <- gwas_dt %>%
		pivot_longer(
			cols = tidyselect::starts_with(log10p_prefix),
			names_to = "ancestry",
			values_to = "log10p",
			values_drop_na = TRUE
		)

	lead_row <- gwas_dt_long[which.max(gwas_dt_long$log10p), ]
	return(lead_row$rsid)
}


# Build a data.table of credible set variants.
# mapping_mode controls how credible-set rsIDs are joined to coordinates/metadata:
# - "gwas": merge SuSiE variants (as snp) to mapping_dt$rsid/position.
# - "pip":  merge SuSiE variants (as snp) to mapping_dt$snp and keep all mapping columns.
# - "none": return only expanded credible-set rows without merging.
build_credible_set_dt <- function(mapping_dt = NULL, cred_set_file,
								  mapping_mode = c("gwas", "pip", "none")) {
	mapping_mode <- match.arg(mapping_mode)

	# SuSiE credible sets can be empty for a locus; return an empty table early.
	cred_dt <- fread(cred_set_file)

	if (nrow(cred_dt) == 0) {
		return(data.table())
	}

	# Keep component as character so labels like "1" are preserved exactly.
	cred_dt <- cred_dt %>%
		mutate(component = as.character(component))

	if ("passed_filter" %in% colnames(cred_dt)) {
		# Accept both logical and string/integer encodings of TRUE.
		cred_dt <- cred_dt %>%
			filter(as.logical(passed_filter) | passed_filter %in% c("TRUE", "True", "true", 1, "1"))
	}

	# Expand comma-delimited variant rsIDs to one row per variant.
	# Do not rely on variant_indices because upstream z files may be cropped/reordered.
	cred_long <- cred_dt %>%
		select(component, variants, num_variants) %>%
		mutate(
			variants = strsplit(as.character(variants), ",")
		) %>%
		mutate(snp = purrr::map(variants, ~ trimws(.x))) %>%
		select(component, num_variants, snp) %>%
		unnest(snp)

	if (identical(mapping_mode, "none")) {
		return(as.data.table(cred_long))
	}

	if (is.null(mapping_dt)) {
		stop("mapping_dt is required when mapping_mode is not 'none'")
	}

	if (identical(mapping_mode, "gwas")) {
		gwas_with_rsid <- copy(as.data.table(mapping_dt))
		if (!all(c("rsid", "position") %in% colnames(gwas_with_rsid))) {
			stop("mapping_dt must contain rsid and position columns for mapping_mode = 'gwas'")
		}
		gwas_with_rsid <- gwas_with_rsid %>%
			transmute(snp = as.character(rsid), position)

		return(
			merge(
				as.data.table(cred_long),
				as.data.table(gwas_with_rsid),
				by = "snp",
				all.x = TRUE,
				allow.cartesian = TRUE
			) %>%
				as.data.table()
		)
	}

	# mapping_mode == "pip"
	pip_with_rsid <- copy(as.data.table(mapping_dt))
	if (!"snp" %in% colnames(pip_with_rsid)) {
		stop("mapping_dt must contain a snp column for mapping_mode = 'pip'")
	}
	pip_with_rsid[, snp := as.character(snp)]

	return(
		merge(
			as.data.table(cred_long),
			pip_with_rsid,
			by = "snp",
			all.x = TRUE,
			allow.cartesian = TRUE
		) %>%
			as.data.table()
	)
}


plot_gwas_tracks <- function(pval_df, chrom, start, end,
							 cred_set_file,
							 x_breaks = NULL, annotate_rsids = NULL,
							 annotate_positions = NULL,
							 ld_data_list = NULL,
							 leading_snp = NULL,
							 mode = c("ivw", "mama"),
							 repel_kws = list(),
							 ...) {
	mode <- match.arg(mode)

	# Restrict GWAS points to the requested genomic window.
	plot_dt <- as.data.table(pval_df)
	plot_dt <- plot_dt[position >= start & position <= end]

	# Harmonized across modes: use GWAS-based merge for credible-set coordinates.
	cred_plot_dt <- build_credible_set_dt(plot_dt, cred_set_file, mapping_mode = "gwas")
	cred_plot_dt <- cred_plot_dt[position >= start & position <= end]

	# Build one track per ancestry for faceted plotting, including LD scores if available.
	eur_dt <- plot_dt %>%
		transmute(
			rsid,
			position,
			ancestry = "EUR",
			log10p = log10p_eur
		)
	if (!is.null(ld_data_list[["EUR"]])) {
		eur_dt <- eur_dt %>% left_join(ld_data_list[["EUR"]], by = "rsid")
	}

	afr_dt <- plot_dt %>%
		transmute(
			rsid,
			position,
			ancestry = "AFR",
			log10p = log10p_afr
		)
	if (!is.null(ld_data_list[["AFR"]])) {
		afr_dt <- afr_dt %>% left_join(ld_data_list[["AFR"]], by = "rsid")
	}

	thai_dt <- plot_dt %>%
		transmute(
			rsid,
			position,
			ancestry = "THAI",
			log10p = log10p_thai
		)
	if (!is.null(ld_data_list[["THAI"]])) {
		thai_dt <- thai_dt %>% left_join(ld_data_list[["THAI"]], by = "rsid")
	}

	track_dt <- bind_rows(eur_dt, afr_dt, thai_dt) %>%
		filter(!is.na(log10p)) %>%
		mutate(
			in_credible_set = FALSE,
			ld_stat = if ("ld_stat" %in% names(.)) ld_stat else NA_real_
		)

	# Mark credible set variants with their component membership.
	if (nrow(cred_plot_dt) > 0) {
		if (!"snp" %in% colnames(cred_plot_dt) && "variant" %in% colnames(cred_plot_dt)) {
			cred_plot_dt <- cred_plot_dt %>% mutate(snp = variant)
		}

		cred_plot_tbl <- tibble::as_tibble(cred_plot_dt) %>%
			transmute(
				rsid = snp,
				position,
				cs_component = paste0("CS", component),
				cs_variant_count = num_variants,
				pip = if ("pip" %in% colnames(cred_plot_dt)) pip else NA_real_
			) %>%
			arrange(desc(pip)) %>%
			distinct(rsid, position, .keep_all = TRUE)

		track_dt <- track_dt %>%
			left_join(cred_plot_tbl, by = c("rsid", "position")) %>%
			mutate(
				in_credible_set = !is.na(cs_component),
				cs_component = if_else(is.na(cs_component), "Background", cs_component)
			)
	} else {
		track_dt <- track_dt %>%
			mutate(
				cs_component = "Background",
				cs_variant_count = NA_integer_,
				pip = NA_real_
			)
	}

	# Label variants on all ancestry panels.
	# Explicitly requested rsIDs are always labeled, even if they are not in a
	# credible set. In addition, all SNPs from small credible sets (<5 variants)
	# are labeled even if not requested.
	label_base_dt <- track_dt %>%
		filter(ancestry == "EUR") %>%
		distinct(
			rsid,
			position,
			ancestry,
			log10p,
			in_credible_set,
			cs_component,
			cs_variant_count,
			.keep_all = TRUE
		)

	annotate_rsids_clean <- character(0)
	if (!is.null(annotate_rsids) && length(annotate_rsids) > 0) {
		annotate_rsids_clean <- unique(trimws(as.character(annotate_rsids)))
		annotate_rsids_clean <- annotate_rsids_clean[nchar(annotate_rsids_clean) > 0]
	}

	label_dt <- label_base_dt %>%
		mutate(is_requested = rsid %in% annotate_rsids_clean) %>%
		filter(
			is_requested |
				(in_credible_set & !is.na(cs_variant_count) & cs_variant_count <= 5)
		) %>%
		mutate(
			label = case_when(
				is_requested ~ rsid,
				!is.na(cs_variant_count) & cs_variant_count < 5 ~ rsid,
				TRUE ~ NA_character_
			),
			ancestry = "EUR"
		) %>%
		filter(!is.na(label))

	p_gwas <- ggplot(track_dt, aes(x = position, y = log10p)) +
		geom_point(
			aes(color = (ld_stat)^2),
			size = 1.5,
			alpha = 0.7
		) +
		# Harmonized across modes: follow regional scripts facet behavior.
		facet_grid(rows = vars(ancestry), scales = "free_y") +
		scale_color_distiller(
			name = bquote("LD (" * r^2 * ") with " * .(leading_snp)),
			palette = "Spectral",
			na.value = "grey70",
			limits = c(0, 1)
		) +
		scale_x_continuous(
			limits = c(start, end),
			breaks = x_breaks,
			labels = function(x) x / 1e6,
			expand = expansion(mult = c(0.01, 0.01))
		) +
		# Keep data coordinates unchanged while preventing marker clipping at panel edges.
		coord_cartesian(clip = "off") +
		labs(x = paste0("chr", chrom, " position"), y = expression(-log[10](p))) +
		theme_minimal() +
		theme(
			panel.grid.minor = element_blank(),
			strip.background = element_rect(fill = "grey95", color = "grey70"),
			legend.position = "bottom"
		)

	# Add hollow circles/squares for credible set variants (different shapes for different components).
	cs_dt <- track_dt %>% filter(in_credible_set)
	if (nrow(cs_dt) > 0) {
		unique_cs <- sort(unique(cs_dt$cs_component[cs_dt$cs_component != "Background"]))
		shape_palette <- c(1, 0, 2, 5, 6, 15)
		shape_map <- setNames(
			shape_palette[seq_len(min(length(unique_cs), length(shape_palette)))],
			unique_cs[seq_len(min(length(unique_cs), length(shape_palette)))]
		)

		p_gwas <- p_gwas +
			geom_point(
				data = cs_dt,
				aes(shape = cs_component),
				fill = NA,
				color = "black",
				size = 1.75,
				stroke = 0.6
			) +
			scale_shape_manual(
				name = "Credible Set",
				values = shape_map,
				guide = "none"
			)
	}

	# Draw vertical guide lines for designated rsIDs across all ancestry panels.
	if (!is.null(annotate_positions) && length(annotate_positions) > 0) {
		p_gwas <- p_gwas +
			geom_vline(
				xintercept = annotate_positions,
				linetype = "dashed",
				color = "#2b8cbe",
				linewidth = 0.35,
				alpha = 0.7
			)
	}

	if (nrow(label_dt) > 0) {
		default_repel <- list(
			data = label_dt,
			mapping = aes(label = label),
			size = 3,
			box.padding = 0.5,
			point.padding = 0.5,
			min.segment.length = 0.2,
			force = 2,
			force_pull = 1,
			show.legend = FALSE,
			color = "black",
			segment.size = 0.3,
			segment.alpha = 0.6
		)
		dots_args <- list(...)
		text_args <- utils::modifyList(default_repel, repel_kws)
		text_args <- utils::modifyList(text_args, dots_args)
		p_gwas <- p_gwas + do.call(geom_text_repel, text_args)
	}

	return(p_gwas)
}


# Extract a display cell-type label from an ATAC track filename (BigWig or bedgraph).
# Handles Corces-2016-style names (e.g. "CD14-3.merge.s20.w150sw.bw") as well as
# generic patterns (e.g. "HSC.ATAC.hg19.bw"); falls back to the bare file stem.
get_cell_type <- function(filename, include_bm_atac = TRUE) {
	bn <- basename(filename)
	if (isTRUE(include_bm_atac)) {
		ct <- sub("-(\\d+|BM-ATAC)\\.merge\\.s20\\.w150sw\\.bw$", "", bn, ignore.case = TRUE)
	} else {
		ct <- sub("-[0-9]+\\.merge\\.s20\\.w150sw\\.bw$", "", bn, ignore.case = TRUE)
	}
	if (identical(ct, bn)) ct <- sub("\\.ATAC\\.hg19\\.(bw|bedgraph)$", "", bn, ignore.case = TRUE)
	if (identical(ct, bn)) ct <- tools::file_path_sans_ext(bn)
	return(ct)
}


# Extract ATAC signal for a specified genomic region from a BigWig or bedgraph file.
extract_atac_region <- function(atac_file, chrom, start, end, include_bm_atac = TRUE) {
	region_start <- start
	region_end <- end
	region_chrom <- normalize_chr(chrom)
	cell_type <- get_cell_type(atac_file, include_bm_atac = include_bm_atac)
	ext <- tolower(tools::file_ext(sub("\\.gz$", "", atac_file, ignore.case = TRUE)))

	if (ext %in% c("bw", "bigwig")) {
		# --- BigWig path ---
		chrom_candidates <- unique(c(as.character(chrom), paste0("chr", region_chrom), region_chrom))
		bw_file <- rtracklayer::BigWigFile(atac_file)
		bw_seqlevels <- tryCatch(names(seqinfo(bw_file)), error = function(e) character())
		query_chrom <- chrom_candidates[1]
		if (length(bw_seqlevels) > 0) {
			matched <- chrom_candidates[chrom_candidates %in% bw_seqlevels]
			if (length(matched) > 0) query_chrom <- matched[1]
		}
		gr <- rtracklayer::import(
			bw_file,
			which = GenomicRanges::GRanges(query_chrom, IRanges::IRanges(region_start, region_end))
		)
		bed_dt <- as.data.table(data.frame(
			chrom = as.character(GenomicRanges::seqnames(gr)),
			start = GenomicRanges::start(gr),
			end = GenomicRanges::end(gr),
			value = gr$score
		))
	} else {
		# --- bedgraph path (plain or .gz) ---
		# bedgraph format: chrom  chromStart  chromEnd  dataValue (0-based half-open)
		raw <- fread(atac_file, header = FALSE, col.names = c("chrom", "start", "end", "value"),
								 colClasses = c("character", "integer", "integer", "numeric"))
		# Convert to 1-based closed for consistency with BigWig handling.
		raw[, start := start + 1L]
		bed_dt <- raw[
			normalize_chr(chrom) == region_chrom &
				end >= region_start &
				start <= region_end
		]
	}
	return(list(cell_type = cell_type, bed_dt = bed_dt))
}


plot_bigwig_track <- function(atac_file, chrom, start, end,
							  x_breaks = NULL, annotate_positions = NULL,
							  include_bm_atac = TRUE,
							  bar_fill = "gray20") {
	region <- extract_atac_region(
		atac_file = atac_file,
		chrom = chrom,
		start = start,
		end = end,
		include_bm_atac = include_bm_atac
	)
	bed_dt <- region$bed_dt
	cell_type <- region$cell_type
	region_start <- start
	region_end <- end
	region_chrom <- normalize_chr(chrom)

	# Clip segments to the requested window and compute bar geometry.
	bed_dt <- as.data.table(bed_dt) %>%
		filter(
			normalize_chr(chrom) == region_chrom,
			end >= region_start,
			start <= region_end
		) %>%
		mutate(
			plot_start = pmax(start, region_start),
			plot_end = pmin(end, region_end),
			midpoint = (plot_start + plot_end) / 2,
			width = pmax(plot_end - plot_start + 1, 1)
		)

	p_bg <- ggplot(bed_dt, aes(x = midpoint, y = value)) +
		geom_col(aes(width = width), fill = bar_fill, color = NA) +
		scale_x_continuous(limits = c(start, end), breaks = x_breaks,
									expand = expansion(mult = c(0.01, 0.01))) +
		labs(y = cell_type) +
		theme_minimal() +
		theme(
			panel.grid.minor = element_blank(),
			axis.line = element_line(color = 'black', linewidth = 0.25),
			plot.title = element_text(size = 10)
		)

	if (!is.null(annotate_positions) && length(annotate_positions) > 0) {
		p_bg <- p_bg +
			geom_vline(xintercept = annotate_positions, linetype = "dashed",
								color = "#2b8cbe", linewidth = 0.35, alpha = 0.7)
	}

	return(p_bg)
}
