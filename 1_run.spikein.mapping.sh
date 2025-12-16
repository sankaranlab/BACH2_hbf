#!/bin/bash

lib=${HOME}/source
software=${HOME}/software
trimmomaticbin=${HOME}/software/Trimmomatic-0.36
adapterpath=${HOME}/source/adapter
bowtie2bin=${HOME}/conda_env/HiC/bin
samtoolsbin=${HOME}/conda_env/HiC/bin
chrom_size=${HOME}/Reference_Genome/hg38/hg38.chrom.sizes
raw=${HOME}/group
ecoli=${HOME}/Reference_Genome/Escherichia_coli_K_12_MG1655/NCBI/2001-10-15/Sequence/Bowtie2Index
path=${HOME}
trimdir2=${HOME}/trimmed3
logdir=${HOME}/logs
aligndir=${HOME}/aligned.aug10
spikealign=${HOME}/spikein
workdir=$path

base=$1


# Map Ecoli

($bowtie2bin/bowtie2 -p 2 --dovetail --phred33 -x $ecoli/genome -1 $trimdir2/"$base"_1.paired.fastq.gz -2 $trimdir2/"$base"_2.paired.fastq.gz) 2> $logdir/"$base".bowtie2 | $samtoolsbin/samtools view -bS - > $spikealign/"$base"_aligned_reads.bam

echo $base `$samtoolsbin/samtools view -c -F 12 $spikealign/"$base"_aligned_reads.bam` >> spikein.read 

