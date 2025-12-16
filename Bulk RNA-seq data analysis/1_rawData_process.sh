#!/bin/bash

# step1: prepare input raw RNA-seq data 

expt=$1
name=$expt
software=${HOME}/software
samtools=$software/samtools
STAR=$software/STAR
fq1=raw_data/${expt}_R1_001.fastq.gz
fq2=raw_data/${expt}_R2_001.fastq.gz
outdir=output/$expt

# Here is index of reference genome which is used to do mapping
index=${HOME}/Reference_Genome/hg38/STAR_index/

# step2 : maping your raw reads

$STAR --runThreadN 8 \
     --runMode alignReads \
     --genomeDir "$index" \
     --readFilesIn "$fq1" "$fq2" \
     --readFilesCommand zcat \
     --outSAMtype BAM SortedByCoordinate \
     --chimSegmentMin 20 \
     --outFileNamePrefix $outdir/$expt

# step3: sort bam file

$samtools view -F 516 -b Aligned.out.sam | samtools sort -@ 12 - > $name.sorted.bam

# Here is gene annotation file (can be download from UCSC)

gtf=/broad/sankaranlab/xww/Reference_Genome/hg38/Homo_sapiens.GRCh38.113.geneName.gtf

# step4: count raw reads for each gene

$software/featureCounts -p -B -a $gtf -t exon -g gene_name -o $outdir/$name.counts.txt $outdir/$name.sorted.bam

cut -f1,7 $outdir/$name.counts.txt > $outdir/$name.count.bed

wait
exit;

