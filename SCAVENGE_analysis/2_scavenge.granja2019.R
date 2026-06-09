#!/broad/sankaranlab/ajlee/conda_envs/sc_env/bin/R


library(Seurat)
library(Signac)
library(SCAVENGE)
library(chromVAR)
library(gchromVAR)
library(SummarizedExperiment)
library(BSgenome.Hsapiens.UCSC.hg19)
library(dplyr)


dir.create('scav_res')
dir.create('seurat_objects')

set.seed(1234)

#
variant_file <- 'input_data/hbf_meta.finemap.hg19.bed'

#
atac_data <- readRDS('/broad/sankaranlab/ajlee/projects/b_all/data/scmultiome_Granja2019/1_collect_data/processed_data/normal_bm.peaks.granja2019.rds')
count_mat <- atac_data@assays$data$counts

peak_data <- data.frame(atac_data@rowRanges)[, c('seqnames','start','end')]
peak_list <- paste0(peak_data$seqnames,':', peak_data$start, '-', peak_data$end)

rownames(count_mat) <- peak_list
atac_obj <- CreateSeuratObject(counts=count_mat, assay = 'peaks', meta.data=data.frame(atac_data@colData))
atac_obj <- RunTFIDF(atac_obj)

##===========================================================================================================================
## scavenge
df <- atac_obj@meta.data[,c('orig.ident','Clusters','BioClassification')]
peak_list <- makeGRangesFromDataFrame(peak_data)

#
se_obj <- SummarizedExperiment(assays=list(counts=atac_obj@assays$peaks@data), rowRanges=peak_list, colData=df)
assayNames(se_obj) <- 'counts'

se_obj <- addGCBias(se_obj, genome=BSgenome.Hsapiens.UCSC.hg19)
peak_bg <- getBackgroundPeaks(object = se_obj, niterations = 200)

## import trait
trait_import <- importBedScore(rowRanges(se_obj), variant_file, colidx=5)

dev <- computeWeightedDeviations(
	object = se_obj,
	weights = trait_import,
	background_peaks = peak_bg
)

## reformat results
zscore_vec <- as.vector(t(SummarizedExperiment::assays(dev)[['z']]))
seed_idx <- seedindex(zscore_vec, percent_cut=0.05)
scale_factor <- cal_scalefactor(z_score=zscore_vec, percent_cut = 0.01)

## construct m-knn graph
peak_by_cell_mat <- SummarizedExperiment::assay(se_obj)
tfidf_mat <- tfidf(
	bmat=peak_by_cell_mat,
	mat_binary=TRUE,
	TF=TRUE,
	log_TF=TRUE
)

lsi_mat <- do_lsi(mat = tfidf_mat, dims = 30)
mutualknn30 <- getmutualknn(lsimat = lsi_mat, num_k = 30)

## network propagation
np_score <- randomWalk_sparse(intM = mutualknn30, queryCells = rownames(mutualknn30)[seed_idx], gamma = 0.05)
omit_idx <- np_score==0
print(sum(omit_idx))

mutualknn30 <- mutualknn30[!omit_idx, !omit_idx]
np_score <- np_score[!omit_idx]

#
trs_raw <- capOutlierQuantile(x = np_score, q_ceiling = 1) |> max_min_scale()
trs_raw <- trs_raw * scale_factor

trs_obj <- subset(atac_obj, cells=names(trs_raw))
trs_obj$trs_raw <- trs_raw

#
trs_cap_01 <- capOutlierQuantile(x = np_score, q_ceiling = 0.999) |> max_min_scale()
trs_cap_01 <- trs_cap_01 * scale_factor
trs_obj$trs_cap_01 <- trs_cap_01

#
trs_cap_1 <- capOutlierQuantile(x = np_score, q_ceiling = 0.99) |> max_min_scale()
trs_cap_1 <- trs_cap_1 * scale_factor
trs_obj$trs_cap_1 <- trs_cap_1

#
trs_cap_5 <- capOutlierQuantile(x = np_score, q_ceiling = 0.95) |> max_min_scale()
trs_cap_5 <- trs_cap_5 * scale_factor
trs_obj$trs_cap_5 <- trs_cap_5

##
trs_data <- data.frame(barcode=rownames(trs_obj@meta.data), trs_raw=trs_obj$trs_raw, trs_cap_01=trs_obj$trs_cap_01, trs_cap_1=trs_obj$trs_cap_1, trs_cap_5=trs_obj$trs_cap_5, celltype=trs_obj$BioClassification)
write.table(trs_data, paste0('scav_res/hbf_meta.finemap.trs.txt'), row.names=F, col.names=T, sep='\t', quote=F)

##===========================================================================================================================



