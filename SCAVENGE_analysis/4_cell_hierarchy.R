#!/broad/sankaranlab/ajlee/conda_envs/sc_env/bin/R


library(Seurat)
library(Signac)
library(ggplot2)
library(RColorBrewer)
library(viridis)


dir.create('plots')


atac_obj <- readRDS('seurat_objects/normal_bm.atac.granja2019.rds')
trs_data <- read.table('scav_res/hbf_meta.finemap.trs.txt', header=T, row.names=1)

#
trs_obj <- subset(atac_obj, cells=rownames(trs_data))
trs_obj <- AddMetaData(trs_obj, metadata=trs_data[,c('trs_raw','trs_cap_01','trs_cap_1','trs_cap_5')])

#
celltype_list <- c(
'01_HSC',
'02_Early.Eryth',
'03_Late.Eryth',
'04_Early.Baso',
'05_CMP.LMPP',
'06_CLP.1',
'07_GMP',
'08_GMP.Neut',
'09_pDC',
'10_cDC', 
'11_CD14.Mono.1',
'12_CD14.Mono.2',
#'13_Unk',
#'14_Unk',
'15_CLP.2',
'16_Pre.B',
'17_B',
'18_Plasma',
'19_CD8.N',
'20_CD4.N1', 
'21_CD4.N2',
'22_CD4.M',
'23_CD8.EM',
'24_CD8.CM',
'25_NK'
#'26_Unk'
)

trs_obj <- subset(trs_obj, subset= BioClassification %in% celltype_list) 

##=============================================================================================
## umap

embedding <- data.frame(trs_obj@meta.data[, c('UMAP1','UMAP2')])
trs_obj$celltype <- factor(trs_obj$BioClassification, levels=celltype_list)

## cell type
coord1 <- tapply(embedding$UMAP1, trs_obj$celltype, mean)
coord2 <- tapply(embedding$UMAP2, trs_obj$celltype, mean)
text_df <- data.frame(coord1, coord2)

celltype_colors <- colorRampPalette(brewer.pal(8, 'Accent'))(length(levels(trs_obj$celltype)))

pdf('plots/atac.granja_2019.celltype.umap.pdf', width=5, height=5, pointsize=3)
plot(
	embedding$UMAP1, embedding$UMAP2,
	col=celltype_colors[trs_obj$celltype],
	xlim=c(min(embedding[,'UMAP1']), max(embedding[,'UMAP1'])),
	ylim=c(min(embedding[,'UMAP2']), max(embedding[,'UMAP2'])),
	xlab='UMAP1', ylab='UMAP2',
	pch=16, cex=0.6, axes=F, ann=F
)
text(x=text_df$coord1, y=text_df$coord2, labels= rownames(text_df), cex=0.6)
dev.off()

png('plots/atac.granja_2019.celltype.umap.png', width=8, height=8, res=1000, units='in')
plot(
	embedding$UMAP1, embedding$UMAP2,
	col=celltype_colors[trs_obj$celltype],
	xlim=c(min(embedding[,'UMAP1']), max(embedding[,'UMAP1'])),
	ylim=c(min(embedding[,'UMAP2']), max(embedding[,'UMAP2'])),
	xlab='UMAP1', ylab='UMAP2',
	pch=16, cex=0.6, axes=F, ann=F
)
dev.off()


## trs
nlength <- 100
score_colors <- viridis(nlength)
col_vec <- score_colors[cut(trs_obj$trs_cap_01, breaks=nlength)]

pdf('plots/atac.granja_2019.trs.umap.pdf', width=5, height=5, pointsize=3)
plot(
	embedding$UMAP1, embedding$UMAP2,
	col=col_vec,
	xlim=c(min(embedding[,'UMAP1']), max(embedding[,'UMAP1'])),
	ylim=c(min(embedding[,'UMAP2']), max(embedding[,'UMAP2'])),
	xlab='UMAP1', ylab='UMAP2',
	pch=16, cex=0.6, axes=F, ann=F
)
dev.off()

png('plots/atac.granja_2019.trs.umap.png', width=8, height=8, res=1000, units='in')
plot(
	embedding$UMAP1, embedding$UMAP2,
	col=col_vec,
	xlim=c(min(embedding[,'UMAP1']), max(embedding[,'UMAP1'])),
	ylim=c(min(embedding[,'UMAP2']), max(embedding[,'UMAP2'])),
	xlab='UMAP1', ylab='UMAP2',
	pch=16, cex=0.6, axes=F, ann=F
)
dev.off()


##=============================================================================================
## boxplot
#
my_theme <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
	theme(panel.background = element_rect(fill='white', color='black', linetype='solid')) +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	theme(legend.position='bottom')

#
trs_data <- trs_data[trs_data$celltype %in% celltype_list, ]
trs_data$celltype <- factor(trs_data$celltype, levels=celltype_list)

g1 <- ggplot(trs_data) +
	geom_boxplot(aes(x=celltype, y=trs_cap_01, fill=celltype), outlier.colour=NA) +
	scale_fill_manual(values=celltype_colors) +
	labs(y='trs', x=NULL) +
	my_theme

pdf(paste0('plots/atac.granja_2019.hbf_meta.trs.box.pdf'), height=4, width=6, pointsize=3)
plot(g1)
dev.off()

##=============================================================================================
##
# Set up coordinates
coords_df <- data.frame(
	cell_label = c("HSC", "MPP", "LMPP", "CLP", "GMP-A", "GMP-B", "GMP-C", "CMP", "MEP", "NK", "CD4", "CD8", "B", "pDC", "Mono", "mDC", "Ery", "Mega"),
	x = c( 0,     0,      -5,    -5,      0,        -2,    2,       5,     7,    -10,   -8,    -6,   -4,  -2,     2,     4,      8,     10), 
	y = c(10,     8,      7,     5,       6,        5,     5,       7,     5,     2,     2,     2,    2,   2,     2,     2,      2,     2)
)

#
mean_vec <- tapply(trs_data$trs_cap_1, trs_data$celltype, mean)





#
#score_colors <- viridis(nlength)
#
#df <- gchromvar
#plotdf <- merge(coords_df, df[df$Trait == plottrait, ], by.x = 'cell_label', by.y = 'celltype')
#plotdf$zscore <- as.numeric(plotdf$zscore)
#plotdf$pvalue <-  pnorm(plotdf$zscore, lower.tail = FALSE)
#
#p1 <- ggplot(plotdf, aes(x = x, y = y, color =  -log10(pvalue))) + 
#	geom_point(size = 11) + pretty_plot() +
#	geom_text(aes(label=cell_label),hjust=0.5, vjust=3) + 
#	scale_color_gradientn(colors = THE_PALETTE, name = "-log10(pvalue)") +
#	scale_y_continuous(limits = c(0, 11)) + ggtitle(plottrait)
#
#ggsave(p1, filename = paste0("../output/gchromVAR/rawPDFs/", plottrait, ".pdf"),
#plot(p1)
#height = 8, width = 10)







##=============================================================================================

