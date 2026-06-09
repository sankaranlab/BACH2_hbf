library(susieR)
library(data.table)
# library(Rfast)

args = commandArgs(trailingOnly=TRUE)

z_filename = args[1]
ld_filename = args[2]
prefix = args[3]
n_cred = as.numeric(args[4])

# spit out usage if arguments are missing
if (length(args) < 4) {
  stop("Usage: Rscript susieR.r <z_matrix_file> <ld_matrix_file> <prefix> <n_cred>", call.=FALSE)
}	

if (interactive()) {
  # get current path
	parent <- getwd()
	# output_path = file.path(parent, "L2_CARMA") <--- change to your own output path
	output_path = file.path(parent, "single_anc_finemap")
	#dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
} else {
	output_path = "single_anc_finemap"
}

if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}

# Z-scores and LD matrix
# ld_matrix = paste(parent, ld_filename,sep='/')
# z_matrix = paste(parent, z_filename,sep='/')
# R <- as.matrix(fread(ld_matrix,sep=' ', header = FALSE))
# z_scores <- fread(z_matrix, sep=' ', header=TRUE)
R <- as.matrix(fread(ld_filename,sep='\t', header = FALSE))
# Make sure the LD matrix is symmetric
R[lower.tri(R)] <- t(R)[lower.tri(R)]

z_scores <- fread(z_filename, sep='\t', header=TRUE)
# get z
z_scores$z <- z_scores$beta / z_scores$se

# Convert beta/se to z if needed
# Run SuSiE
res <- susie_rss(z=z_scores$z, R=R, n=median(z_scores$n), L = n_cred)

# Collect result
cred_sets <- unique(unlist(res$sets$cs))
print(length(res$sets$cs))
print(names(res$sets$cs))
str(res$sets$cs)

pip <- res$pip
z_scores$pip=pip
z_scores$cred=0L

if(length(cred_sets)!=0){
	cs_names <- names(res$sets$cs)
	print(cs_names)
	for(L in seq_along(res$sets$cs)){
		cred_id <- L
		if(!is.null(cs_names) && !is.na(cs_names[L])){
			name_id <- suppressWarnings(as.integer(sub('^L', '', cs_names[L])))
			if(!is.na(name_id)) cred_id <- name_id
		}
		cs_idx <- as.integer(res$sets$cs[[L]])
		cs_idx <- cs_idx[cs_idx >= 1 & cs_idx <= nrow(z_scores)]
		if(length(cs_idx) > 0){
			z_scores[cs_idx, cred := as.integer(cred_id)]
		}
	}
	creds=z_scores[z_scores$cred!=0L,]
}else{
	creds <- z_scores[0]
}
out_path <- file.path(output_path, paste0(prefix, ".out"))
fwrite(z_scores, out_path, sep = "\t")
# only write credible set if it's non-empty
if(nrow(creds) > 0){
	out_path_cred <- file.path(output_path, paste0(prefix, ".cred"))
	fwrite(creds, out_path_cred, sep = "\t")
}
