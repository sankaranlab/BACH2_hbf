#!/bin/bash

# Output/Error file (%x=job-name, %j=job-id)
#SBATCH --output=logs/%x.o%j
#SBATCH --error=logs/%x.e%j

# Partition name
#SBATCH --partition=bch-compute

#SBATCH --time=2:00:00
#SBATCH --mem=16G

# Number of nodes
#SBATCH --nodes=1

# Number of tasks (processes)
#SBATCH --ntasks=1               

# Number of CPU cores per task
#SBATCH --cpus-per-task=4

# Email notifications (set to your own address)
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

set -euo pipefail

# Set working directory.
# NOTE: the paths below are environment-specific; adapt them to your system.
LABSHARE="/path/to/shared"
WDIR="${LABSHARE}/BACH2_hbf/StatGen_analysis/MultiSuSiE"
cd $WDIR

THAI_GWAS="${WDIR}/qc-ed/Thai_chrALL_MAF0p1pct_hg38_info6.clean.tsv.gz"
CHAIN="${LABSHARE}/ref_genomes/human/hg38ToHg19.over.chain.gz"

# Extract CHR and BP from gwas, form a bed for liftover
echo "Preparing GWAS summary stats for liftover"
zcat ${THAI_GWAS} | awk 'BEGIN{OFS="\t"} NR>1{print "chr"$2, $3-1, $3}' > ${WDIR}/raw/hg19_sumstats/Thai_gwas_hg38_info6.bed

# liftover
echo ""
echo "$(date) Performing liftover"
liftOver -minMatch=0.98 -multiple ${WDIR}/raw/hg19_sumstats/Thai_gwas_hg38_info6.bed $CHAIN ${WDIR}/raw/hg19_sumstats/Thai_gwas_hg19_info6.bed ${WDIR}/raw/hg19_sumstats/Thai_gwas_info6_unmapped.bed

# merge them and verify:
echo ""
echo "$(date) Merging liftover results"
python scripts/1_merge_liftover_coord.py \
  --gwas-file ${THAI_GWAS} -from hg38 --to hg19 \
  --liftover-bed ${WDIR}/raw/hg19_sumstats/Thai_gwas_hg19_info6.bed \
  -o ${WDIR}/raw/hg19_sumstats/Thai_gwas_hg19_info6_cleaned_merged.tsv.gz \
  --filtered-file ${WDIR}/raw/hg19_sumstats/Thai_gwas_hg19_info6_dbsnp-filtered.tsv.gz \
  --parallel $SLURM_CPUS_PER_TASK

echo ""
echo "Pipeline finished $(date)"
echo "==================================================================="
echo ""

# report job stats
sacct -j $SLURM_JOB_ID --format=JobID,JobName,Partition,NodeList,Start,End,Elapsed,State,ExitCode,MaxRSS,MaxVMSize,AveCPU,AveRSS,AveVMSize
