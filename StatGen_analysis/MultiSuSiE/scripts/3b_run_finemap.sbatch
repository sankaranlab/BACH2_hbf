#!/bin/bash

#SBATCH --output=logs/%x.o%j
#SBATCH --error=logs/%x.e%j
#SBATCH --time=6:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1

# notify me when done (set to your own address)
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=<your-email>

# load module
module load finemap

## take $GENE from sbatch
GENE="${GENE:?GENE is not set. Pass with: sbatch --export=GENE=<GENE> ...}"
: "${POP:=}"

# Set working directory.
# NOTE: the paths below are environment-specific; adapt them to your system.
LABSHARE="/path/to/shared"
WDIR="${LABSHARE}/projects/xhcheng/HbF/bach2_multi_finemap"
cd "$WDIR"

## commands here
z_input=$WDIR/loci/finemap_input/${GENE}_flank1Mb_info6.single.${POP}.z
ld_input=$WDIR/loci/finemap_input/${GENE}_flank1Mb_info6.single.${POP}.ld
prefix=${GENE}_flank1Mb_info6_susie

master_file=$WDIR/loci/finemap_input/${GENE}_flank1Mb_info6.single.${POP}.master

finemap --sss \
  --in-files $WDIR/$master_file \
  --n-causal-snps 20


# report job stats
sacct -j $SLURM_JOB_ID --format=JobID,JobName,Partition,NodeList,Start,End,Elapsed,State,ExitCode,MaxRSS,MaxVMSize,AveCPU,AveRSS,AveVMSize
