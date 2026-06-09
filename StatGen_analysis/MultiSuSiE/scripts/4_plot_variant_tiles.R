library(tidyverse)
library(ggpubr)
library(data.table)
library(ggrepel)
library(RColorBrewer)
library(rtracklayer)
library(patchwork)

b37_GENE_REF <- '/lab-share/Hem-Sankaran-e2/Public/ref_genomes/human/hg19_ncbiRefSeq_RefSeqAll_2026UCSC.tsv.gz'

ATAC_DIR_WENG2024   <- '/lab-share/Hem-Sankaran-e2/Public/inhouse_datasets/bone_marrow_atac_peak_bigwig/hg19_liftOver/'
ATAC_DIR_CORCES2016 <- '/lab-share/Hem-Sankaran-e2/Public/inhouse_datasets/Corces_2016_hg19bw/'

# Load shared step4 utility functions from the sibling scripts utility folder.
get_script_dir <- function() {
    file_arg <- sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))])
    if (length(file_arg) > 0) {
        return(dirname(normalizePath(file_arg[1])))
    }
    return(getwd())
}
source(file.path(get_script_dir(), "..", "2052490848", "step4_utility.R"))


resolve_variant_atac <- function(variant_dt, atac_dt_list) {
    # variant_dt: data.table with columns chrom, position, rsid
    # atac_dt_list: list of lists with cell_type and bed_dt slots
    variant_dt <- copy(variant_dt)
    variant_dt[, chrom := normalize_chr(chrom)]
    variant_dt[, `:=`(start = as.integer(position), end = as.integer(position))]
    setkey(variant_dt, chrom, start, end)

    for (atac_info in atac_dt_list) {
        cell_type <- atac_info$cell_type
        source_name <- atac_info$source_name
        bed_dt <- copy(atac_info$bed_dt)
        atac_col <- make_atac_column_name(source_name, cell_type)
        if (nrow(bed_dt) == 0) {
            variant_dt[, (atac_col) := NA_real_]
            next
        }
        bed_dt[, chrom := normalize_chr(chrom)]
        setkey(bed_dt, chrom, start, end)

        variant_dt[, (atac_col) := {
            hits <- foverlaps(
                .SD[, .(chrom, start, end, rsid)],
                bed_dt,
                by.x = c("chrom", "start", "end"),
                by.y = c("chrom", "start", "end"),
                nomatch = 0L
            )
            if (nrow(hits) > 0) {
                # If multiple hits (e.g. overlapping peaks), take the max signal.
                max_signal <- tapply(hits$value, hits$rsid, max, na.rm = TRUE)
                return(max_signal[rsid])
            } else {
                return(NA_real_)
            }
        }]
    }
    variant_dt[, c("start", "end") := NULL]
    return(variant_dt)
}


parse_args <- function(args) {
    opts <- list()
    i <- 1L
    while (i <= length(args)) {
        key <- args[[i]]
        if (!startsWith(key, "--")) {
            stop(sprintf("Unexpected argument: %s", key))
        }
        key <- sub("^--", "", key)

        if (grepl("=", key, fixed = TRUE)) {
            kv <- strsplit(key, "=", fixed = TRUE)[[1]]
            opts[[kv[1]]] <- kv[2]
            i <- i + 1L
            next
        }

        if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
            opts[[key]] <- TRUE
            i <- i + 1L
        } else {
            if (identical(key, "cell-types")) {
                j <- i + 1L
                vals <- character()
                while (j <= length(args) && !startsWith(args[[j]], "--")) {
                    vals <- c(vals, args[[j]])
                    j <- j + 1L
                }
                opts[[key]] <- paste(vals, collapse = " ")
                i <- j
            } else {
                opts[[key]] <- args[[i + 1L]]
                i <- i + 2L
            }
        }
    }
    return(opts)
}


parse_atac_sources <- function(atac_dir_arg) {
    if (identical(atac_dir_arg, "Weng2024")) {
        return(data.table(source = "Weng2024", dir = ATAC_DIR_WENG2024))
    }
    if (identical(atac_dir_arg, "Corces2016")) {
        return(data.table(source = "Corces2016", dir = ATAC_DIR_CORCES2016))
    }
    if (identical(atac_dir_arg, "both")) {
        return(data.table(
            source = c("Weng2024", "Corces2016"),
            dir = c(ATAC_DIR_WENG2024, ATAC_DIR_CORCES2016)
        ))
    }
    if (!dir.exists(atac_dir_arg)) {
        stop(sprintf("ATAC directory not found: %s", atac_dir_arg))
    }
    return(data.table(source = "Custom", dir = atac_dir_arg))
}


parse_cell_types <- function(cell_types_arg) {
    if (is.null(cell_types_arg) || !nzchar(trimws(cell_types_arg))) {
        return(NULL)
    }

    arg_str <- trimws(cell_types_arg)
    known_sources <- c("Weng2024", "Corces2016", "Custom")
    parsed <- list()

    # Find all source positions in the string
    all_positions <- numeric()
    all_sources <- character()

    for (src_name in known_sources) {
        matches <- gregexpr(paste0(src_name, "\\s*:"), arg_str, ignore.case = TRUE)[[1]]
        if (matches[1] != -1) {
            all_positions <- c(all_positions, matches)
            all_sources <- c(all_sources, rep(src_name, length(matches)))
        }
    }

    if (length(all_positions) == 0) {
        return(NULL)
    }

    # Sort by position
    order_idx <- order(all_positions)
    all_positions <- all_positions[order_idx]
    all_sources <- all_sources[order_idx]

    # For each source, extract everything after its colon until the next source or end
    for (i in seq_along(all_positions)) {
        src_name <- all_sources[i]
        src_pos <- all_positions[i]

        # Find where the colon ends (start of types)
        colon_match <- regexpr(paste0(src_name, "\\s*:"), substr(arg_str, src_pos, nchar(arg_str)), ignore.case = TRUE)
        types_start <- src_pos + attr(colon_match, "match.length")

        # Find end: next source position or end of string
        types_end <- nchar(arg_str) + 1
        if (i < length(all_positions)) {
            types_end <- all_positions[i + 1L]
        }

        # Extract and parse types
        types_str <- trimws(substr(arg_str, types_start, types_end - 1))
        types_str <- trimws(gsub('^"|"$', "", types_str))

        if (!nzchar(types_str)) {
            next
        }

        types <- trimws(strsplit(types_str, ",", fixed = TRUE)[[1]])
        types <- types[nzchar(types)]

        if (length(types) > 0) {
            parsed[[length(parsed) + 1L]] <- data.table(source = src_name, cell_type = types)
        }
    }

    if (length(parsed) == 0) {
        return(NULL)
    }

    out <- rbindlist(parsed, use.names = TRUE)
    out[, source := as.character(source)]
    out[, cell_type := as.character(cell_type)]
    out <- out[nzchar(trimws(source)) & nzchar(trimws(cell_type))]
    out <- unique(out)

    if (nrow(out) == 0) {
        return(NULL)
    }
    return(out)
}



make_atac_column_name <- function(source_name, cell_type) {
    paste0("ATAC_", source_name, "__", cell_type)
}


get_atac_source_fill_scale <- function(source_name) {
    if (identical(source_name, "Weng2024")) {
        return(scale_fill_gradient(low = "#f0f0f0", high = "#b07db7", na.value = "#f0f0f0", name = "ATAC\nWeng2024"))
    }
    if (identical(source_name, "Corces2016")) {
        return(scale_fill_gradient(low = "#f0f0f0", high = "#4f9a94", na.value = "#f0f0f0", name = "ATAC\nCorces2016"))
    }
    return(scale_fill_gradient(low = "#f0f0f0", high = "#b07db7", na.value = "#f0f0f0", name = "ATAC"))
}


parse_pip_sources <- function(pip_arg) {
    pip_arg <- trimws(as.character(pip_arg))
    if (!nzchar(pip_arg)) {
        stop("--pip is empty")
    }

    if (!grepl(":", pip_arg, fixed = TRUE)) {
        return(data.table(file = pip_arg, label = "PIP"))
    }

    parts <- strsplit(pip_arg, ",", fixed = TRUE)[[1]]
    parts <- trimws(parts)
    parts <- parts[nzchar(parts)]
    if (length(parts) == 0) {
        stop("--pip mapping string is empty")
    }

    parsed <- lapply(parts, function(part) {
        # Try to match with single quotes first
        m <- regexec("^(.+?):\\s*'([^']+)'\\s*$", part)
        g <- regmatches(part, m)[[1]]
        
        # If no match, try without quotes
        if (length(g) == 0) {
            m <- regexec("^(.+?):\\s*([^\\s]+)\\s*$", part)
            g <- regmatches(part, m)[[1]]
        }
        
        if (length(g) != 3) {
            stop(sprintf("Invalid --pip entry: %s. Expected format <file>:'label' or <file>:label", part))
        }
        list(file = trimws(g[2]), label = trimws(g[3]))
    })

    out <- rbindlist(lapply(parsed, as.data.table), use.names = TRUE)
    if (anyDuplicated(out$label)) {
        stop("Duplicate labels in --pip mapping are not allowed")
    }
    return(out)
}


read_and_merge_pip_sources <- function(pip_sources, credible_snps_dt, chrom_arg = NULL) {
    seed_dt <- unique(as.data.table(credible_snps_dt)[, .(rsid = as.character(rsid))])
    if (nrow(seed_dt) == 0) {
        stop("No credible-set SNPs to annotate with PIP values")
    }
    seed_dt[, `:=`(position = as.integer(NA), chrom = NA_character_)]

    pip_cols <- character()
    pip_label_map <- character()
    safe_label_names <- make.unique(make.names(pip_sources$label, allow_ = TRUE))

    for (i in seq_len(nrow(pip_sources))) {
        src_file <- pip_sources$file[[i]]
        src_label <- pip_sources$label[[i]]
        if (!file.exists(src_file)) {
            stop(sprintf("pip file not found: %s", src_file))
        }

        dt <- as.data.table(fread(src_file))
        if (!("rsid" %in% names(dt)) && ("snp" %in% names(dt))) {
            setnames(dt, "snp", "rsid")
        }
        if (!("rsid" %in% names(dt))) stop(sprintf("pip file missing rsid/snp: %s", src_file))
        if (!("position" %in% names(dt))) stop(sprintf("pip file missing position: %s", src_file))
        if (!("pip" %in% names(dt))) stop(sprintf("pip file missing pip: %s", src_file))

        dt <- dt[as.character(rsid) %chin% seed_dt$rsid]

        chrom_col <- NULL
        for (nm in c("chrom", "chr", "chromosome")) {
            if (nm %in% names(dt)) {
                chrom_col <- nm
                break
            }
        }

        pip_col <- paste0("PIP_", safe_label_names[[i]])
        pip_label_map[[pip_col]] <- src_label
        pip_cols <- c(pip_cols, pip_col)

        if (nrow(dt) == 0) {
            seed_dt[, (pip_col) := NA_real_]
            next
        }

        if (is.null(chrom_col)) {
            dt <- dt[, .(
                rsid = as.character(rsid),
                position_src = as.integer(position),
                chrom_src = NA_character_,
                pip_src = as.numeric(pip)
            )]
        } else {
            dt <- dt[, .(
                rsid = as.character(rsid),
                position_src = as.integer(position),
                chrom_src = normalize_chr(get(chrom_col)),
                pip_src = as.numeric(pip)
            )]
        }

        # If a file has duplicate rsids, retain the row with the largest PIP for that source.
        dt <- dt[order(-pip_src)]
        dt <- dt[, .SD[1], by = rsid]
        setnames(dt, "pip_src", pip_col)

        seed_dt <- merge(seed_dt, dt, by = "rsid", all.x = TRUE, sort = FALSE)
        seed_dt[, position := data.table::fcoalesce(position, position_src)]
        seed_dt[, chrom := data.table::fcoalesce(chrom, chrom_src)]
        seed_dt[, c("position_src", "chrom_src") := NULL]
    }

    if (!is.null(chrom_arg)) {
        seed_dt[is.na(chrom) | chrom == "", chrom := normalize_chr(chrom_arg)]
    }
    if (all(is.na(seed_dt$chrom) | seed_dt$chrom == "")) {
        stop("No chromosome column found in any pip file. Provide --chrom.")
    }
    if (all(is.na(seed_dt$position))) {
        stop("No positions found for credible-set SNPs in the supplied pip file(s).")
    }

    return(list(pip_dt = seed_dt, pip_cols = pip_cols, pip_label_map = pip_label_map))
}


usage <- function() {
    msg <- c(
        "Usage:",
        "Rscript step4_plot_variant_tiles.R \\",
        "  --credible-set path/to/*_credible_sets.tsv \\",
        "  --pip \"file1.tsv:'label1',file2.tsv:'label2'\" \\",
        "  --atac-dir Weng2024|Corces2016|both|/path/to/atac_tracks \\",
        "  --out output.png \\",
        "  [--chrom 6] [--atac-pattern regex] [--cell-types Weng2024:\"type1,type2,...,typen\" Corces2016:\"type1,type2,...,typen\"] [--rank-var-by position|pip|atac] [--split-pip] [--title text] [--width 6] [--height 12] [--dpi 300]",
        "",
        "Notes:",
        "- The script plots variants in passed_filter == TRUE credible sets.",
        "- If pip file lacks a chromosome column, --chrom is required.",
        "- --pip supports either one file path or a comma-separated mapping string: file:'label'.",
        "- Each mapped PIP source is shown as its own PIP tile column.",
        "- --atac-dir both plots Weng2024 and Corces2016 as separate ATAC panels with separate color scales.",
        "- Each ATAC track/cell type is shown as its own tile column.",
        "- --cell-types specifies cell types per source (space-separated source:\"type1,type2,...\" entries).",
        "  Example: --cell-types Weng2024:\"HSC,MPP\" Corces2016:\"CD34,CD14\"",
        "- PIP is colored with a pseudo-log scale.",
        "- split-pip isolates PIP == 1 into its own top color.",
        "- rank-var-by controls row order from top to bottom."
    )
    cat(paste(msg, collapse = "\n"), "\n")
}


find_atac_files <- function(atac_dir, atac_pattern = NULL) {
    files <- list.files(atac_dir, full.names = TRUE)
    files <- files[grepl("\\.(bw|bigwig|bedgraph)(\\.gz)?$", files, ignore.case = TRUE)]
    if (!is.null(atac_pattern) && !identical(atac_pattern, "")) {
        files <- files[grepl(atac_pattern, basename(files), ignore.case = TRUE)]
    }
    return(sort(files))
}


row_max_na <- function(x) {
    if (is.null(dim(x))) {
        return(as.numeric(x))
    }
    out <- apply(x, 1, function(v) {
        if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)
    })
    return(as.numeric(out))
}


compute_plot_dimensions <- function(n_rows, n_atac_cols, n_pip_cols,
                                    width_arg = NULL, height_arg = NULL) {
    n_cols <- n_atac_cols + n_pip_cols
    if (n_rows <= 0 || n_cols <= 0) {
        return(list(width = 6, height = 6))
    }

    # Tile geometry and non-panel framing space (axis labels, title margin, legends).
    tile_size_in <- 0.16
    extra_w <- 3.2
    extra_h <- 1.6

    if (is.null(width_arg) && is.null(height_arg)) {
        width <- max(6, n_cols * tile_size_in + extra_w)
        height <- max(6, n_rows * tile_size_in + extra_h)
        return(list(width = width, height = height))
    }

    if (!is.null(width_arg) && is.null(height_arg)) {
        body_w <- max(width_arg - extra_w, 1)
        body_h <- body_w * (n_rows / n_cols)
        height <- max(6, body_h + extra_h)
        return(list(width = width_arg, height = height))
    }

    if (is.null(width_arg) && !is.null(height_arg)) {
        body_h <- max(height_arg - extra_h, 1)
        body_w <- body_h * (n_cols / n_rows)
        width <- max(6, body_w + extra_w)
        return(list(width = width, height = height_arg))
    }

    return(list(width = width_arg, height = height_arg))
}


plot_single_atac_panel <- function(plot_dt, atac_cols, source_name, show_variant_labels = TRUE) {
    atac_tile <- data.table::melt(
        copy(plot_dt),
        id.vars = intersect(c("rsid", "position", "rsid_factor", "atac_rank_sum"), names(plot_dt)),
        measure.vars = atac_cols,
        variable.name = "metric",
        value.name = "value"
    )
    atac_tile[, metric := sub(paste0("^ATAC_", source_name, "__"), "", metric)]
    atac_tile[, metric := factor(metric, levels = sub(paste0("^ATAC_", source_name, "__"), "", atac_cols))]

    base_theme <- theme_bw(base_size = 10) +
        theme(
            panel.grid = element_blank(),
            axis.title = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
        )

    p <- ggplot(atac_tile, aes(x = metric, y = rsid_factor, fill = value)) +
        geom_tile(color = "black", linewidth = 0.25) +
        coord_fixed(ratio = 1) +
        scale_x_discrete(position = "top", expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0)) +
        get_atac_source_fill_scale(source_name) +
        base_theme +
        theme(
            axis.text.x = element_blank(),
            axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0),
            axis.ticks.x = element_blank(),
            axis.ticks.x.top = element_blank()
        )

    if (!show_variant_labels) {
        p <- p +
            theme(
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank()
            )
    }

    return(p)
}


make_pip_fill_scale <- function(pip_values, split_pip = FALSE) {
    pip_values <- as.numeric(pip_values)
    pip_values <- pip_values[!is.na(pip_values)]
    split_pip_color <- "#203a8f"

    if (length(pip_values) == 0) {
        return(
            scale_fill_gradient(
                low = "#f0f0f0",
                high = "#cc1f5a",
                na.value = "#f0f0f0",
                name = "log(PIP)",
                # trans = scales::pseudo_log_trans(base = 10)
            )
        )
    }

    if (!split_pip) {
        return(
            scale_fill_gradient(
                low = "#f0f0f0",
                high = "#cc1f5a",
                na.value = "#f0f0f0",
                name = "log(PIP)",
                # trans = scales::pseudo_log_trans(base = 10)
            )
        )
    }

    non_one <- pip_values[pip_values < 1]
    if (length(non_one) == 0) {
        return(
            scale_fill_gradient(
                low = "#f7d6df",
                high = split_pip_color,
                na.value = "#f0f0f0",
                name = "log(PIP)"
            )
        )
    }

    min_non_one <- min(non_one, na.rm = TRUE)
    max_non_one <- max(non_one, na.rm = TRUE)

    if (!is.finite(min_non_one) || !is.finite(max_non_one) || max_non_one <= 0) {
        return(
            scale_fill_gradient(
                low = "#f0f0f0",
                high = "#cc1f5a",
                na.value = "#f0f0f0",
                name = "PIP",
                trans = scales::pseudo_log_trans(base = 10)
            )
        )
    }

    split_top <- max_non_one * 1.15
    if (split_top <= max_non_one) {
        split_top <- max_non_one + max(max_non_one * 0.15, 1e-6)
    }

    return(
        scale_fill_gradientn(
            colours = c("#f0f0f0", "#cc1f5a", split_pip_color),
            values = c(0, 0.88, 1),
            limits = c(min_non_one, split_top),
            na.value = "#f0f0f0",
            name = "PIP",
            trans = scales::pseudo_log_trans(base = 10),
            breaks = c(min_non_one, max_non_one, split_top),
            labels = c(signif(min_non_one, 2), signif(max_non_one, 2), "1")
        )
    )
}


plot_variant_tiles <- function(plot_dt, title = NULL, rank_var_by = "position", split_pip = FALSE,
                               pip_label_map = NULL, atac_source_order = NULL,
                               selected_cell_types = NULL) {
    plot_dt <- copy(plot_dt)
    atac_cols <- grep("^ATAC_", names(plot_dt), value = TRUE)
    pip_cols <- grep("^PIP_", names(plot_dt), value = TRUE)
    if (length(atac_cols) == 0) {
        stop("plot_dt must contain one or more ATAC_* columns")
    }
    if (length(pip_cols) == 0) {
        if ("pip" %in% names(plot_dt)) {
            plot_dt[, PIP_PIP := pip]
            pip_cols <- "PIP_PIP"
        } else {
            stop("plot_dt must contain one or more PIP_* columns")
        }
    }

    atac_col_info <- data.table(
        col = atac_cols,
        source = sub("^ATAC_([^_]+)__.*$", "\\1", atac_cols),
        cell_type = sub("^ATAC_[^_]+__", "", atac_cols)
    )

    # Enforce cell-type filtering at plot time so only nominated cell types are shown.
    if (!is.null(selected_cell_types) && nrow(selected_cell_types) > 0) {
        keep_idx <- rep(FALSE, nrow(atac_col_info))
        for (src in unique(atac_col_info$source)) {
            src_types <- selected_cell_types[tolower(source) == tolower(src), cell_type]
            if (length(src_types) == 0) {
                next
            }
            keep_idx <- keep_idx | (atac_col_info$source == src & atac_col_info$cell_type %in% src_types)
        }
        atac_col_info <- atac_col_info[keep_idx]
        atac_cols <- atac_col_info$col
        if (length(atac_cols) == 0) {
            stop("No ATAC columns remain after --cell-types filtering.")
        }
    }
    if (is.null(atac_source_order)) {
        atac_source_order <- unique(atac_col_info$source)
    }

    valid_rankings <- c("position", "pip", "atac")
    if (!(rank_var_by %in% valid_rankings)) {
        stop(sprintf("rank_var_by must be one of: %s", paste(valid_rankings, collapse = ", ")))
    }

    if (rank_var_by == "position") {
        if (!("position" %in% names(plot_dt))) {
            stop("plot_dt must contain position when rank_var_by = 'position'")
        }
        plot_dt <- plot_dt[order(position, na.last = TRUE)]
    } else if (rank_var_by == "pip") {
        pip_mat <- as.matrix(plot_dt[, ..pip_cols])
        plot_dt[, pip_rank_value := apply(pip_mat, 1, function(x) {
            if (all(is.na(x))) {
                return(NA_real_)
            }
            return(max(x, na.rm = TRUE))
        })]
        plot_dt <- plot_dt[order(-pip_rank_value, na.last = TRUE)]
    } else if (rank_var_by == "atac") {
        atac_mat <- as.matrix(plot_dt[, ..atac_cols])
        plot_dt[, atac_rank_sum := apply(atac_mat, 1, function(x) {
            if (all(is.na(x))) {
                return(NA_real_)
            }
            return(sum(x, na.rm = TRUE))
        })]
        plot_dt <- plot_dt[order(-atac_rank_sum, na.last = TRUE)]
    }

    plot_dt[, rsid := as.character(rsid)]
    plot_dt[, rsid_factor := factor(rsid, levels = rev(rsid))]
    pip_long <- data.table::melt(
        copy(plot_dt),
        id.vars = intersect(c("rsid", "position", "rsid_factor", "atac_rank_sum", "pip_rank_value"), names(plot_dt)),
        measure.vars = pip_cols,
        variable.name = "metric",
        value.name = "pip_raw"
    )
    pip_long[, pip_plot := pip_raw]

    if (split_pip) {
        non_one <- pip_long[!is.na(pip_raw) & pip_raw < 1, pip_raw]
        if (length(non_one) > 0) {
            max_non_one <- max(non_one, na.rm = TRUE)
            split_top <- max_non_one * 1.15
            if (split_top <= max_non_one) {
                split_top <- max_non_one + max(max_non_one * 0.15, 1e-6)
            }
            pip_long[pip_raw >= 1, pip_plot := split_top]
        }
    }

    pip_tile <- copy(pip_long)
    if (!is.null(pip_label_map) && length(pip_label_map) > 0) {
        pip_metric_labels <- vapply(pip_cols, function(x) {
            if (!is.null(pip_label_map[[x]])) pip_label_map[[x]] else sub("^PIP_", "", x)
        }, character(1))
    } else {
        pip_metric_labels <- sub("^PIP_", "", pip_cols)
    }
    pip_tile[, metric := factor(metric, levels = pip_cols, labels = pip_metric_labels)]

    base_theme <- theme_bw(base_size = 10) +
        theme(
            panel.grid = element_blank(),
            axis.title = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
        )

    p_pip <- ggplot(pip_tile, aes(x = metric, y = rsid_factor, fill = log10(pip_plot))) +
        geom_tile(color = "black", linewidth = 0.25) +
        coord_fixed(ratio = 1) +
        scale_x_discrete(position = "top", expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0)) +
        make_pip_fill_scale(pip_tile$pip_plot, split_pip = split_pip) +
        base_theme +
        theme(
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.x = element_blank(),
            axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0),
            axis.ticks.x = element_blank(),
            axis.ticks.x.top = element_blank()
        )
    # Overplot PIP==1 tiles with a fixed colour so the split is guaranteed to render
    # correctly regardless of how scale_fill_gradientn resolves values at its upper limit.
    if (split_pip) {
        pip_one_data <- pip_tile[!is.na(pip_raw) & pip_raw >= 1]
        if (nrow(pip_one_data) > 0) {
            p_pip <- p_pip +
                geom_tile(data = pip_one_data, fill = "#203a8f", color = "black", linewidth = 0.25)
        }
    }

    atac_panels <- list()
    width_vec <- numeric()

    for (source_idx in seq_along(atac_source_order)) {
        source_name <- atac_source_order[[source_idx]]
        source_cols <- atac_col_info[source == source_name, col]
        if (length(source_cols) == 0) next

        show_variant_labels <- !(length(atac_source_order) > 1 && source_idx > 1)
        atac_panels[[source_name]] <- plot_single_atac_panel(
            plot_dt,
            source_cols,
            source_name,
            show_variant_labels = show_variant_labels
        )
        width_vec <- c(width_vec, length(source_cols))
    }

    bottom_panels <- c(atac_panels, list(PIP = p_pip))
    bottom_widths <- c(width_vec, length(pip_cols))

    p <- patchwork::wrap_plots(bottom_panels, widths = bottom_widths, guides = "collect") &
        theme(legend.position = "right")

    if (!is.null(title) && !identical(title, "")) {
        p <- p + patchwork::plot_annotation(title = title)
    }
    return(p)
}



if(! interactive()) {
    message("[step4] Starting variant tile plot generation...")
    args <- parse_args(commandArgs(trailingOnly = TRUE))

    if (isTRUE(args$help) || isTRUE(args$h)) {
        usage()
        quit(save = "no", status = 0)
    }
    message("[step4] Parsed arguments successfully")
    print(args)

    cred_set_file <- args[["credible-set"]]
    if (is.null(cred_set_file)) cred_set_file <- args[["cred-set"]]
    pip_file <- args[["pip"]]
    atac_dir <- args[["atac-dir"]]
    out_file <- args[["out"]]

    required_missing <- c(
        "credible-set" = is.null(cred_set_file),
        "pip" = is.null(pip_file),
        "atac-dir" = is.null(atac_dir),
        "out" = is.null(out_file)
    )
    if (any(required_missing)) {
        missing_names <- names(required_missing)[required_missing]
        stop(sprintf("Missing required args: %s", paste(missing_names, collapse = ", ")))
    }

    atac_sources <- parse_atac_sources(atac_dir)
    message(sprintf("[step4] Found %d ATAC source(s): %s", nrow(atac_sources), paste(atac_sources$source, collapse = ", ")))
    
    selected_cell_types <- parse_cell_types(args[["cell-types"]])
    if (!is.null(selected_cell_types)) {
        message(sprintf("[step4] Filtering to cell types: %s", paste(unique(selected_cell_types$cell_type), collapse = ", ")))
    } else {
        message("[step4] No cell-type filter applied")
    }

    if (!file.exists(cred_set_file)) stop(sprintf("credible set file not found: %s", cred_set_file))

    message("[step4] Loading credible set from:", cred_set_file)
    cred_long <- build_credible_set_dt(cred_set_file = cred_set_file, mapping_mode = "none")
    if (nrow(cred_long) == 0) {
        message("[step4] No variants in passed-filter credible sets")
        p_empty <- ggplot() + theme_void() + ggtitle("No variants in passed-filter credible sets")
        ggsave(out_file, p_empty, width = 6, height = 4, dpi = 300, bg = "white")
        message(sprintf("[step4] Wrote empty plot: %s", out_file))
        quit(save = "no", status = 0)
    }
    message(sprintf("[step4] Loaded %d credible-set entries", nrow(cred_long)))

    credible_snps <- unique(cred_long[, .(rsid = as.character(snp))])
    message("[step4] Loading PIP scores from:", pip_file)
    pip_sources <- parse_pip_sources(pip_file)
    pip_loaded <- read_and_merge_pip_sources(
        pip_sources = pip_sources,
        credible_snps_dt = credible_snps,
        chrom_arg = args[["chrom"]]
    )
    pip_dt <- pip_loaded$pip_dt
    pip_cols <- pip_loaded$pip_cols
    pip_label_map <- pip_loaded$pip_label_map
    message(sprintf("[step4] Loaded PIP for %d variants, %d columns", nrow(pip_dt), length(pip_cols)))

    message("[step4] Merging credible sets with PIP scores and genomic positions...")
    variant_dt <- merge(
        cred_long[, .(component, num_variants, rsid = as.character(snp))],
        pip_dt[, c("rsid", "chrom", "position", pip_cols), with = FALSE],
        by = "rsid",
        all.x = TRUE,
        allow.cartesian = TRUE
    )

    variant_dt <- variant_dt[!is.na(position)]
    if (nrow(variant_dt) == 0) stop("No variants with valid genomic positions after merge.")
    message(sprintf("[step4] Merged to %d variants at valid positions", nrow(variant_dt)))

    # Union of TRUE credible-set variants; if repeated across components, keep max PIP entry.
    variant_dt[, pip_rank_value := apply(as.matrix(.SD), 1, function(x) {
        if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
    }), .SDcols = pip_cols]
    variant_dt <- variant_dt[order(-pip_rank_value)]
    variant_dt <- variant_dt[, .SD[1], by = rsid]

    region_chr <- unique(variant_dt$chrom)
    if (length(region_chr) != 1) {
        stop("Variants span multiple chromosomes; this script currently expects one locus/chromosome.")
    }
    region_start <- min(variant_dt$position, na.rm = TRUE)
    region_end   <- max(variant_dt$position, na.rm = TRUE)
    message(sprintf("[step4] Plotting region: %s:%d-%d", region_chr, region_start, region_end))

    message("[step4] Loading ATAC tracks...")
    atac_dt_list <- list()
    for (i in seq_len(nrow(atac_sources))) {
        source_name <- atac_sources$source[[i]]
        source_dir <- atac_sources$dir[[i]]
        source_files <- find_atac_files(source_dir, args[["atac-pattern"]])
        if (length(source_files) == 0) {
            next
        }

        # Filter by cell types if specified
        if (!is.null(selected_cell_types)) {
            source_cell_types <- vapply(source_files, get_cell_type, character(1))
            
            # Get cell types for this source
            source_spec_types <- selected_cell_types[tolower(source) == tolower(source_name), cell_type]
            
            source_files <- source_files[source_cell_types %in% source_spec_types]
            if (length(source_files) == 0) {
                next
            }
        }

        message(sprintf("Reading %d %s ATAC tracks for %s:%d-%d", length(source_files), source_name, region_chr, region_start, region_end))
        source_dt_list <- lapply(source_files, function(f) {
            atac_info <- extract_atac_region(f, chrom = region_chr, start = region_start, end = region_end)
            atac_info$source_name <- source_name
            atac_info
        })
        atac_dt_list <- c(atac_dt_list, source_dt_list)
    }
    if (length(atac_dt_list) == 0) {
        stop("No ATAC files found after applying source and tree filters.")
    }
    message(sprintf("[step4] Loaded %d ATAC track(s)", length(atac_dt_list)))

    variant_atac <- resolve_variant_atac(
        variant_dt = variant_dt[, .(rsid, chrom, position)],
        atac_dt_list = atac_dt_list
    )
    atac_cols <- grep("^ATAC_", names(variant_atac), value = TRUE)
    if (length(atac_cols) == 0) stop("No ATAC columns were created from the supplied tracks.")

    # Reorder ATAC columns by cell type if specified
    if (!is.null(selected_cell_types)) {
        atac_col_info <- data.table(
            col = atac_cols,
            source = sub("^ATAC_([^_]+)__.*$", "\\1", atac_cols),
            cell_type = sub("^ATAC_[^_]+__", "", atac_cols)
        )
        
        # Create order mapping for each source
        atac_col_info[, order := NA_integer_]
        
        for (src in unique(atac_col_info$source)) {
            # Get cell types specified for this source
            src_types <- selected_cell_types[tolower(source) == tolower(src), cell_type]
            src_order <- setNames(seq_along(src_types), src_types)
            atac_col_info[source == src, order := src_order[cell_type]]
        }
        
        # Filter to only specified cell types and reorder
        atac_col_info <- atac_col_info[!is.na(order)]
        setorder(atac_col_info, source, order)
        atac_cols <- atac_col_info$col
    }

    atac_source_order <- atac_sources$source

    plot_dt <- merge(
        variant_atac[, c("rsid", "position", atac_cols), with = FALSE],
        variant_dt[, c("rsid", "position", pip_cols), with = FALSE],
        by = c("rsid", "position"),
        all.x = TRUE
    )
    rank_var_by <- "position"
    split_pip <- isTRUE(args[["split-pip"]])
    title <- args[["title"]]
    width_in <- if (!is.null(args[["width"]])) as.numeric(args[["width"]]) else NULL
    height_in <- if (!is.null(args[["height"]])) as.numeric(args[["height"]]) else NULL
    
    message(sprintf("[step4] Generating plot with %d variants, %d ATAC columns, %d PIP columns", nrow(plot_dt), length(atac_cols), length(pip_cols)))
    dims <- compute_plot_dimensions(
        n_rows = nrow(plot_dt),
        n_atac_cols = length(atac_cols),
        n_pip_cols = length(pip_cols),
        width_arg = width_in,
        height_arg = height_in
    )
    width <- dims$width
    height <- dims$height
    dpi <- if (!is.null(args[["dpi"]])) as.integer(args[["dpi"]]) else 300L

    p <- plot_variant_tiles(
        plot_dt,
        title = title,
        rank_var_by = rank_var_by,
        split_pip = split_pip,
        pip_label_map = pip_label_map,
        atac_source_order = atac_source_order,
        selected_cell_types = selected_cell_types
    )

    message("[step4] Saving plots...")
    for (ext in c("png", "pdf")) { # , "svg"
        out_file <- paste0(tools::file_path_sans_ext(out_file), ".", ext)
        if (ext == "png"){
            ggsave(out_file, p, width = width, height = height, dpi = dpi)
        } else {
            ggsave(out_file, p, width = width, height = height)
        }
        message(sprintf("[step4] Wrote: %s", out_file))
    }
    message("[step4] Complete!")
}