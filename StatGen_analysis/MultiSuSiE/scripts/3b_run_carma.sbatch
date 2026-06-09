#!/bin/bash

#SBATCH --output=logs/%x.o%j
#SBATCH --error=logs/%x.e%j
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1

# notify me when done (set to your own address)
#SBATCH --mail-type=END
#SBATCH --mail-user=<your-email>

# load module
module load r/4.4.2

# Set working directory.
# NOTE: the paths below are environment-specific; adapt them to your system.
LABSHARE="/path/to/shared"
WDIR="${LABSHARE}/projects/xhcheng/HbF/bach2_multi_finemap"
cd "$WDIR"

## take $GENE_LABEL from sbatch
: "${GENE_LABEL:=}"
# : "${POP:=}"

# if gene label is for BACH2, then subPop should be ""/"_subEAS" in the input, "_subCSA"/"_subEAS" in the output; if gene label is for BCL11A or HBS1L, then subPop should be "_subCSA" or "_subEAS"
# if [ "$GENE_LABEL" == "BACH2_flank300kb" ]; then
#     INPUT_SUBPOPS=("" "_subEAS")
#     OUTPUT_SUBPOPS=("_subCSA" "_subEAS")
# else
#     INPUT_SUBPOPS=("_subCSA" "_subEAS")
#     OUTPUT_SUBPOPS=("_subCSA" "_subEAS")
# fi


for POP in EUR AFR; do

        echo "[$(date +'%Y-%m-%d %H:%M:%S')] Running CARMA for ${GENE_LABEL} in ${POP} sub-population ${output_subPop}"
        ## commands here
        z_input=$WDIR/loci/finemap_input/${GENE_LABEL}_info6_mac40.single.${POP}.z
        ld_input=$WDIR/loci/finemap_input/${GENE_LABEL}_info6_mac40.single.${POP}.ld
        prefix=${GENE_LABEL}_info6_mac40_carma.${POP}

        Rscript $WDIR/scripts/carma.r $z_input $ld_input $prefix

        echo ""
done

echo ""
POP="THAI"
# for i in "${!INPUT_SUBPOPS[@]}"; do
    # input_subPop="${INPUT_SUBPOPS[$i]}"
    # output_subPop="${OUTPUT_SUBPOPS[$i]}"
for subPop in {"_subCSA","_subEAS"}; do
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] Running CARMA for ${GENE_LABEL} in ${POP} sub-population ${subPop}"
    ## commands here
    z_input=$WDIR/loci/finemap_input/${GENE_LABEL}_info6_mac40${subPop}.single.${POP}.z
    ld_input=$WDIR/loci/finemap_input/${GENE_LABEL}_info6_mac40${subPop}.single.${POP}.ld
    prefix=${GENE_LABEL}_info6_mac40${subPop}_carma.${POP}

    Rscript $WDIR/scripts/carma.r $z_input $ld_input $prefix

    echo ""
done


# report job stats
sacct -j $SLURM_JOB_ID --format=JobID,JobName,Partition,NodeList,Start,End,Elapsed,State,ExitCode,MaxRSS,MaxVMSize,AveCPU,AveRSS,AveVMSize
