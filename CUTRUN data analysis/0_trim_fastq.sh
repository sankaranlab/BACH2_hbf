lib=${HOME}/source
software=${HOME}/software
trimmomaticbin=${HOME}/software/Trimmomatic-0.36
adapterpath=${HOME}/source/adapter
bowtie2bin=${HOME}/conda_env/HiC/bin
samtoolsbin=${HOME}/conda_env/HiC/bin
chrom_size=${HOME}/Reference_Genome/hg38/hg38.chrom.sizes
raw=${HOME}/group/fastq
#len=`cat length`
len=150
workdir=${HOME}/
trimdir=$workdir/trimmed
trimdir2=$workdir/trimmed3
logdir=$workdir/logs
bt2idx=${HOME}/Reference_Genome/hg38/bowtie_index
aligndir=$workdir/aligned.aug10

#mkdir $trimdir
#mkdir $trimdir2
#mkdir $logdir
#mkdir $aligndir

base=$1

${HOME}/conda_env/HiC/bin/java -jar $trimmomaticbin/trimmomatic-0.36.jar PE -threads 1 -phred33 $raw/"$base"_R1_001.fastq.gz $raw/"$base"_R2_001.fastq.gz $trimdir/"$base"_1.paired.fastq.gz $trimdir/"$base"_1.unpaired.fastq.gz $trimdir/"$base"_2.paired.fastq.gz $trimdir/"$base"_2.unpaired.fastq.gz ILLUMINACLIP:$adapterpath/Truseq3.PE.fa:2:15:4:4:true LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:25

#>&2 echo "Second stage trimming $base ..."
#>&2 date
${HOME}/conda_env/sc_env/bin/python $lib/kseq_python.py $trimdir/"$base"_1.paired.fastq.gz $len $trimdir2/"$base"_1.paired.fastq.gz
${HOME}/conda_env/sc_env/bin/python $lib/kseq_python.py $trimdir/"$base"_2.paired.fastq.gz $len $trimdir2/"$base"_2.paired.fastq.gz

