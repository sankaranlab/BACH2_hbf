#!/bin/bash

#SBATCH --output=logs/%x.o%j
#SBATCH --error=logs/%x.e%j
#SBATCH --time=2:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=1

# notify me when done (set to your own address)
#SBATCH --mail-type=END
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

## load module
module load r/4.4.2

# Set working directory.
# NOTE: the paths below are environment-specific; adapt them to your system.
LABSHARE="/path/to/shared"
WDIR="${LABSHARE}/projects/xhcheng/HbF/bach2_multi_finemap"
cd "$WDIR"

## take $GENE from sbatch
: "${GENE:=}"
: "${FLANK:=300kb}"
# : "${POP:=}"
# : "${subPop:=}"

for POP in {"EUR","AFR"}; do
    echo "$(date) Running SuSiE for ${GENE} in ${POP}"
    ## commands here
    z_input=$WDIR/loci/finemap_input/${GENE}_flank${FLANK}_info6_mac40.single.${POP}.z
    ld_input=$WDIR/loci/finemap_input/${GENE}_flank${FLANK}_info6_mac40.single.${POP}.ld
    prefix=${GENE}_flank${FLANK}_info6_mac40.${POP}.susie
    # take n_cred from sbatch command line argument; if not provided, set to default value 5
    # n_cred=${n_cred:-15}

    # Rscript adapted_scripts/susieR.r $z_input $ld_input $prefix $n_cred
    # Rscript $WDIR/scripts/susie_wDiagnostics.r $z_input $ld_input $prefix $n_cred
    for n_cred in {15,10,5,3}; do
        Rscript $WDIR/scripts/susieR.r $z_input $ld_input "${prefix}_L${n_cred}" $n_cred
    done
    echo ""
done

echo ""
POP="THAI"
for subPop in {"_subCSA","_subEAS"}; do
        echo "$(date) Running SuSiE for ${GENE} in ${POP} sub-population ${subPop}"
    ## commands here
    z_input=$WDIR/loci/finemap_input/${GENE}_flank${FLANK}_info6_mac40${subPop}.single.${POP}.z
    ld_input=$WDIR/loci/finemap_input/${GENE}_flank${FLANK}_info6_mac40${subPop}.single.${POP}.ld
    prefix=${GENE}_flank${FLANK}_info6_mac40${subPop}.${POP}.susie
    # take n_cred from sbatch command line argument; if not provided, set to default value 5
    # n_cred=${n_cred:-15}

    # Rscript adapted_scripts/susieR.r $z_input $ld_input $prefix $n_cred
    # Rscript $WDIR/scripts/susie_wDiagnostics.r $z_input $ld_input $prefix $n_cred
    for n_cred in {15,10,5,3}; do
        Rscript $WDIR/scripts/susieR.r $z_input $ld_input "${prefix}_L${n_cred}" $n_cred
    done
    echo ""
done



# report job stats
sacct -j $SLURM_JOB_ID --format=JobID,JobName,Partition,NodeList,Start,End,Elapsed,State,ExitCode,MaxRSS,MaxVMSize,AveCPU,AveRSS,AveVMSize
