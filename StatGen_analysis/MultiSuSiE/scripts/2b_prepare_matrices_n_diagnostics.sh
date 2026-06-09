#!/bin/bash

# Output/Error file (%x=job-name, %j=job-id)
#SBATCH --output=logs/%x.o%j
#SBATCH --error=logs/%x.e%j

# Partition name
#SBATCH --partition=bch-compute

#SBATCH --time=12:00:00
#SBATCH --mem=32G

# Number of nodes
#SBATCH --nodes=1

# Number of tasks (processes)
#SBATCH --ntasks=1               

# Number of CPU cores per task
#SBATCH --cpus-per-task=1

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
WDIR="${LABSHARE}/projects/xhcheng/HbF/bach2_multi_finemap"
cd $WDIR

# Parameters
GENE=${GENE:-BACH2}
# FLANK="1Mb"
FLANK=${FLANK:-300kb}
ANC=${ANC:-3} # number of ancestries (2 or 3)
DIAG=${DIAG:-1} # whether to generate diagnostics plots
subPOP=${subPOP:-EAS} # sub-population label for Thai LD (e.g. EAS or CSA)
THAI_GWAS="raw/hg19_sumstats/Thai_gwas_hg19_info6_cleaned_merged.tsv.gz"
# THAI_GWAS="raw/og_Thai_chrALL_INFO4MAF0p1pct_hg19.tsv"
EUR_GWAS="fema/hg19_output/METAL_hbf_inv_EUR_info6_mac40_SE_gcOff_hg19.tsv.gz"
AFR_GWAS="fema/hg19_output/METAL_hbf_inv_AFR_info6_mac40_SE_gcOff_hg19.tsv.gz"
OUTPUT_PREFIX="loci/finemap_input/${GENE}_flank${FLANK}_info6_mac40_${ANC}ANC_sub${subPOP}"
LD_DIR="ld/ref_ld/"
# THAI_LD="ld/thai_dbsnp_QCed_${GENE}_1Mb_flanking.ld"
# THAI_LD="ld/ref_ld/UKBB.CSA.ldadj.${GENE}_flank1Mb.npy"
THAI_LD="ld/ref_ld/UKBB.${subPOP}.ldadj.${GENE}_flank1Mb.npy"
FINEMAP_OUTDIR="single_anc_finemap/finemap_output/"

SKIP_DIAG_FLAG=()
if [[ "$DIAG" -eq 0 ]]; then
    SKIP_DIAG_FLAG=(--skip-global-diagnostics)
fi

# Run step2b
if [[ $ANC -eq 3 ]]; then
    python scripts/step2b_prepare_3way_matrices_n_diagnostics.py \
    --gene "$GENE" \
    --flank "$FLANK" \
    --input-files \
        "$THAI_GWAS",THAI \
        "$EUR_GWAS",EUR \
        "$AFR_GWAS",AFR \
    --output-prefix "$OUTPUT_PREFIX" \
    --local-ld-dir "$LD_DIR" \
    --thai-ld-path "$THAI_LD" \
    --finemap-outdir "$FINEMAP_OUTDIR" \
    "${SKIP_DIAG_FLAG[@]}" 
elif [[ $ANC -eq 2 ]]; then
    python scripts/step2b_prepare_2way_matrices_n_diagnostics.py \
    --gene "$GENE" \
    --flank "$FLANK" \
    --input-files \
        "$EUR_GWAS",EUR \
        "$AFR_GWAS",AFR \
    --output-prefix "$OUTPUT_PREFIX" \
    --local-ld-dir "$LD_DIR" \
    "${SKIP_DIAG_FLAG[@]}"
else
    echo "Error: Invalid number of ancestries specified. Use --anc 2 or --anc 3."
    exit 1
fi

# report job stats
sacct -j $SLURM_JOB_ID --format=JobID,JobName,Partition,NodeList,Start,End,Elapsed,State,ExitCode,MaxRSS,MaxVMSize,AveCPU,AveRSS,AveVMSize
