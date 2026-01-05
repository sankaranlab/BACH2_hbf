#!/bin/bash

lib=${HOME}/source

#peak_file=$1 #a narrowPeak file
base=$1
#mbase=`basename $peak_file _peaks.narrowPeak`
memebin=${HOME}/meme/bin
bedopsbin=${HOME}/software/bin
bedtoolsbin=${HOME}/conda_env/scenicplus/bin
genome_sequence=${HOME}/Reference_Genome/hg38/hg38.fa

p=0.05
#motif_dir=${HOME}/sources/JASPER
motif_dir=${HOME}/motif_calling/motif_for_use
#base=`basename $peak_file .narrowPeak`
workdir=${HOME}
#dir=`dirname $peak_file`
fa_dir=${HOME}/motif_calling/filtered.fa
peakfile=${HOME}/sorted/ecoli_norm_remove_dup/peak/"$base".downsampled_peaks.narrowPeak


$bedtoolsbin/bedtools getfasta -fi $genome_sequence -bed $peakfile -fo $fa_dir/"$base".fa
${HOME}/conda_env/sc_env/bin/python $lib/fix_sequence.py $fa_dir/"$base".fa

outdir=${HOME}/motif_calling
#mkdir $outdir/$base
for d in $outdir $outdir/$base; do
if [ ! -d $d ]; then
mkdir $d
fi
done

for m in `ls -1 $motif_dir`; do
motif=`basename $m .meme`
fimo_d=$outdir/$base/fimo2.$motif
if [ ! -d $fimo_d ]; then
mkdir $fimo_d
fi
$memebin/fimo --thresh $p --parse-genomic-coord -oc $fimo_d $motif_dir/"$motif".meme $fa_dir/"$base".fa
$bedopsbin/gff2bed --max-mem 8G < $fimo_d/fimo.gff | awk 'BEGIN {IFS="\t"; OFS="\t";} {print $1,$2,$3,$4,$5,$6}' > $fimo_d/fimo.bed
done

