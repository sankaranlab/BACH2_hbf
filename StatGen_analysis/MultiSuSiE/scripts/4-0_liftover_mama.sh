#!/bin/bash

# Output/Error file (%x=job-name, %j=job-id)
#SBATCH --output=logs/%x.o%j
#SBATCH --error=logs/%x.e%j

# Partition name
#SBATCH --partition=bch-compute

#SBATCH --time=4:00:00
#SBATCH --mem=32G

# Number of nodes
#SBATCH --nodes=1

# Number of tasks (processes)
#SBATCH --ntasks=1               

# Number of CPU cores per task
#SBATCH --cpus-per-task=4

# Email notifications
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=<your-email>

echo "========================= Job Information ========================"
echo "Job ID: $SLURM_JOB_ID"
echo "Job Name: $SLURM_JOB_NAME"
echo "Node: $SLURM_NODELIST"
echo "Start Time: $(date)"
echo "Working Directory: $(pwd)"
echo "Number of CPUs: $SLURM_CPUS_PER_TASK"
echo "Clock Time Limit: $SLURM_TIMELIMIT"
echo "Memory per Node: $SLURM_MEM_PER_NODE"
echo "==================================================================="
echo ""

module load bcftools

set -euo pipefail

# set working directory
LABSHARE="/path/to/shared"
WDIR="${LABSHARE}/BACH2_hbf/StatGen_analysis/MultiSuSiE"
cd $WDIR

CHAIN="${LABSHARE}/ref_genomes/human/hg38ToHg19.over.chain.gz"

# Get $MAMA_GWAS  and $POPfrom sbatch arguments
# Extract CHR and BP from MAMA, form a bed for liftover
# MAMA file header: SNP	CHR	BP	A1	A2	FREQ	BETA	SE	Z	P	N_EFF	N_ORIG
echo "Preparing MAMA summary stats for liftover"
cat ${MAMA_GWAS} | awk 'BEGIN{OFS="\t"} NR>1{print "chr"$2, $3-1, $3}' > ${WDIR}/raw/hg19_sumstats/MAMA_${POP}_hg38.bed

# liftover
echo ""
echo "$(date) Performing liftover"
liftOver -minMatch=0.98 -multiple ${WDIR}/raw/hg19_sumstats/MAMA_${POP}_hg38.bed $CHAIN ${WDIR}/raw/hg19_sumstats/MAMA_${POP}_hg19.bed ${WDIR}/raw/hg19_sumstats/MAMA_${POP}_unmapped.bed

# merge them and verify:
echo ""
echo "$(date) Merging liftover results"
python scripts/step1c_merge_liftover_coord.py \
  --gwas-file ${MAMA_GWAS} -from hg38 --to hg19 \
  --liftover-bed ${WDIR}/raw/hg19_sumstats/MAMA_${POP}_hg19.bed \
  -o ${WDIR}/raw/hg19_sumstats/MAMA_${POP}_hg19_cleaned_merged.tsv.gz \
  --filtered-file ${WDIR}/raw/hg19_sumstats/MAMA_${POP}_hg19_dbsnp-filtered.tsv.gz \
  --parallel $SLURM_CPUS_PER_TASK --debug

echo ""
echo "Pipeline finished $(date)"
echo "==================================================================="
echo ""

# report job stats
sacct -j $SLURM_JOB_ID --format=JobID,JobName,Partition,NodeList,Start,End,Elapsed,State,ExitCode,MaxRSS,MaxVMSize,AveCPU,AveRSS,AveVMSize