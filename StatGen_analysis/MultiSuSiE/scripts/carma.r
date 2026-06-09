library(CARMA)
library(data.table)

args = commandArgs(trailingOnly=TRUE)
# prefix = args[1]
z_filename = args[1]
ld_filename = args[2]
prefix = args[3]

# spit out usage if arguments are missing
if (length(args) < 3) {
  stop("Usage: Rscript carma.r <z_matrix_file> <ld_matrix_file> <prefix>", call.=FALSE)
}	

if (interactive()) {
  # get current path
	parent <- getwd()
	# output_path = file.path(parent, "L2_CARMA") <--- change to your own output path
	output_path = file.path(parent, "single_anc_finemap")
	#dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
	ld_filename = paste(parent, ld_filename,sep='/')
	z_filename = paste(parent, z_filename,sep='/')
} else {
	output_path = "single_anc_finemap"
}


#foreach(i = filename, .packages = c("CARMA", "data.table")) %dopar% {
R <- as.matrix(fread(ld_filename, header = FALSE, sep="\t"))
# Make sure the LD matrix is symmetric
R[lower.tri(R)] <- t(R)[lower.tri(R)]

z_scores <- fread(z_filename, sep='\t', header=TRUE)
z <- z_scores$beta / z_scores$se 
# z <- z_scores$z
# z.list <- list(z)
# ld.list <- list(R)

z.list <- list()
z.list[[1]] <- z
ld.list <- list()
ld.list[[1]] <- as.matrix(R)
# lambda.list <- list(1) # <-- this is just a parameter for CARMA
lambda.list <- list()
lambda.list[[1]] <- 1


result <- CARMA(z.list, ld.list, 
				lambda.list = lambda.list, 
				outlier.switch = TRUE, 
				num.causal = 5, 
				printing.log=TRUE)
pip <- result[[1]]$PIPs
z_scores$pip <- pip
cs = result[[1]]$'Credible set'[[2]]
if(length(cs)!=0){
	creds=z_scores[unlist(cs),]
	creds$cs=1
}else{
	temp = z_scores[order(-z_scores$pip),]
	temp = temp[temp$pip>0.001,]
	temp$cum_pip <- cumsum(temp$pip)
	creds <- temp[temp$cum_pip <= 0.95,]
	creds$cs=0
}

out_path <- file.path(output_path, paste0(prefix, ".out"))
fwrite(z_scores, out_path, sep = "\t")
out_path_cred <- file.path(output_path, paste0(prefix, ".cred"))
fwrite(creds, out_path_cred, sep = "\t")
#}
