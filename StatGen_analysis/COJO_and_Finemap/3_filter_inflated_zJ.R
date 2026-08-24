library(data.table)

args <- commandArgs(trailingOnly = TRUE)
cojo_file <- args[1]
foldchange_cutoff <- as.numeric(args[2])
out_prefix <- args[3]

cat("Using fold-change cutoff of", foldchange_cutoff, "\n")

# Read in the data. Header: Chr     SNP     bp      refA    freq    b       se      p       n       freq_geno       bJ      bJ_se       pJ      LD_r
dt <- fread(cojo_file)

# compte Z and ZJ from b/se and bJ/bJ_se, respectively
dt[, Z  := b / se]
dt[, ZJ := bJ / bJ_se]

# define % increase threshold (e.g., 30% increase)
pct_thresh <- foldchange_cutoff

# identify and print rows that will be removed
dt_removed <- dt[abs(ZJ / Z) > (1 + pct_thresh)]
cat("Rows removed (not meeting threshold):", nrow(dt_removed), "\n")
# print(dt_removed[, .(Chr, SNP, bp, Z, p, ZJ, pJ, LD_r)])

# filter: ZJ is no more than pct_thresh larger than Z
dt_fc <- dt[abs(ZJ / Z) <= (1 + pct_thresh)]
cat("Rows retained (meeting threshold):", nrow(dt_fc), "\n")


# reorder columns for clarity
dt_fc <- dt_fc[, .(Chr, Position=bp, rsID=SNP, refA, freq, beta=b, se, Z, p, 
                    N_eff=n, freq_eff=freq_geno, beta_J=bJ, 
                    se_J=bJ_se, Z_J=ZJ, p_J=pJ, LD_r)]

# write out if needed
ST4_name <- paste0(out_prefix, "_filtered_cutoff", 1+pct_thresh, ".tsv")
cat("Writing filtered results to", ST4_name, "\n")
fwrite(dt_fc, ST4_name, sep = "\t")

