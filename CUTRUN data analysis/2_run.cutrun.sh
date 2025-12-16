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
bt2idx=${HOME}/Reference_Genome/hg38/bowtie_index
path=${HOME}
trimdir2=${HOME}/trimmed3
logdir=${HOME}/logs
aligndir=${HOME}/aligned.aug10
spikealign=${HOME}/spikein
workdir=$path

base=$1

# Step1: map raw read to hg38 reference genome by bowtie2

($bowtie2bin/bowtie2 -p 10 --dovetail --phred33 -x $bt2idx/hg38 -1 $trimdir2/"$base"_1.paired.fastq.gz -2 $trimdir2/"$base"_2.paired.fastq.gz) 2> $logdir/"$base".bowtie2 | $samtoolsbin/samtools view -bS - > $aligndir/"$base"_aligned_reads.bam

$samtoolsbin/samtools sort $path/aligned.aug10/"$base"_aligned_reads.bam  -o $path/sorted/"$base".bam

$samtoolsbin/samtools index $path/sorted/"$base".bam

# Step2: mark and remove duplication

${HOME}/conda_env/HiC/bin/java -Xmx50g -jar ${HOME}/conda_env/HiC/share/picard-2.27.5-0/picard.jar MarkDuplicates \
INPUT=$path/sorted/$base.bam OUTPUT=$path/dedup/$base.bam \
METRICS_FILE=metrics.$base.txt \
REMOVE_DUPLICATES=true
$samtoolsbin/samtools index $path/dedup/$base.bam
#bamCoverage --bam $path/dedup/$base.bam  --binSize 1 --normalizeUsing CPM --outFileName $path/dedup/${base}_all_CPM.bw
# Map Ecoli

# Step3: Peak calling 

outdir=$path/peak
${HOME}/conda_env/scenicplus/bin/macs2 callpeak -t $workdir/dedup/"$base".bam -g hs -f BAMPE -n $base --outdir $outdir -q 0.01 -B --SPMR --keep-dup all 2> $logdir/"$base_file".macs2

# Step4: create monoclonal read

$samtoolsbin/samtools sort -n $path/dedup/$base.bam | $samtoolsbin/samtools view -bf 2  | ${HOME}/conda_env/scenicplus/bin/bedtools bamtobed -bedpe -i stdin | awk -v OFS='\t' '{print $1,$2,$6,".",1}' | perl $lib/sorted_bed_merge_redundant_lines_V2.pl | sort -k1,1 -k2,2n -k3,3n > $workdir/monoclonal/$base.monoclonal.bed

# Step5: calculate genome coverage

${HOME}/conda_env/scenicplus/bin/bedtools genomecov -bg -i $workdir/monoclonal/$base.monoclonal.bed -g $chrom_size > $workdir/monoclonal/$base.monoclonal.bed.bg

# Step6: create bigwig track

$software/bedGraphToBigWig $workdir/monoclonal/$base.monoclonal.bed.bg $chrom_size $workdir/bigwig/$base.bw


