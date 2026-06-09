#!/bin/bash

# Output/Error file (%x=job-name, %j=job-id)
#SBATCH --output=logs/%x.o%j
#SBATCH --error=logs/%x.e%j

# Partition name
#SBATCH --partition=bch-compute

#SBATCH --time=12:00:00
#SBATCH --mem=16G

# Number of nodes
#SBATCH --nodes=1

# Number of tasks (processes)
#SBATCH --ntasks=1               

# Number of CPU cores per task
#SBATCH --cpus-per-task=4


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

THREADS="${SLURM_CPUS_PER_TASK:-1}"
if [[ "$THREADS" -lt 1 ]]; then
    THREADS=1
fi
export LC_ALL=C

# load tools
module load bcftools

# set working directory
LABSHARE="/lab-share/Hem-Sankaran-e2/Public"
WDIR="${LABSHARE}/projects/xhcheng/HbF/bach2_multi_finemap"
cd $WDIR

DBSNP_VCF="/lab-share/Hem-Sankaran-e2/Public/public_datasets/dbSNP/build_157/VCF/GCF_000001405.25.gz"
OUTDIR="${WDIR}/fema/hg19_output"
mkdir -p "$OUTDIR"

if [[ ! -f "$DBSNP_VCF" ]]; then
    echo "ERROR: dbSNP VCF not found: $DBSNP_VCF" >&2
    exit 1
fi

if [[ ! -f "${DBSNP_VCF}.tbi" && ! -f "${DBSNP_VCF}.csi" ]]; then
    echo "ERROR: dbSNP index not found: ${DBSNP_VCF}.tbi or ${DBSNP_VCF}.csi" >&2
    exit 1
fi

# receives filelist from --export
INPUT="${INPUT:-}"
OUTNAME="${OUTNAME:-}"
echo "Received input: ${INPUT}"
echo "Output prefix: ${OUTNAME}"
echo "Using threads: ${THREADS}"

if [[ -z "$INPUT" || -z "$OUTNAME" ]]; then
    echo "ERROR: INPUT and OUTNAME must be provided via --export" >&2
    exit 1
fi

if [[ ! -f "$INPUT" ]]; then
    echo "ERROR: Input file not found: $INPUT" >&2
    exit 1
fi

OUTFILE="${OUTDIR}/${OUTNAME}"

TMPDIR="/temp_work/ch268464/bach2_finemap/step0d_metal_to_hg19_${SLURM_JOB_ID}"
mkdir -p "$TMPDIR"
# trap 'rm -rf "$TMPDIR"' EXIT

read_input() {
    if [[ "$INPUT" == *.gz ]]; then
        gzip -cd "$INPUT"
    else
        cat "$INPUT"
    fi
}

echo "[$(date)] Extracting rsIDs from input"
read_input | awk 'BEGIN{FS=OFS="\t"} NR==1{next} $1!=""{print $1}' | sort --parallel="$THREADS" -u > "$TMPDIR/rsids.txt"

echo "[$(date)] Querying dbSNP with bcftools"
bcftools view --threads "$THREADS" -i "ID=@$TMPDIR/rsids.txt" "$DBSNP_VCF" -Ou \
    | bcftools query -f '%ID\t%CHROM\t%POS\n' \
    | awk 'BEGIN{FS=OFS="\t"}
        # Keep only canonical primary chromosomes for each rsID (1-22, X, Y, MT).
        function normalize_chr(c, n) {
            gsub(/^chr/,"",c)
            if (c ~ /^[0-9]+$/ || c=="X" || c=="Y") return c
            if (c=="M" || c=="MT") return "MT"
            if (c ~ /^NC_/) {
                n = substr(c,4,6) + 0
                if (n>=1 && n<=22) return n ""
                if (n==23) return "X"
                if (n==24) return "Y"
                if (n==12920) return "MT"
            }
            return c
        }
        function chr_rank(c) {
            if (c ~ /^[0-9]+$/) return c + 0
            if (c=="X") return 23
            if (c=="Y") return 24
            if (c=="MT") return 25
            return 999
        }
        {
            id=$1
            chr=normalize_chr($2)
            pos=$3
            r=chr_rank(chr)
            # Ignore all non-primary contigs entirely.
            if (r > 25) next
            # For primary chromosomes, keep lowest rank/position if duplicates exist.
            if (!(id in best_rank) || r < best_rank[id] || (r == best_rank[id] && pos+0 < best_pos[id])) {
                best_rank[id]=r
                best_chr[id]=chr
                best_pos[id]=pos+0
            }
        }
        END {
            for (id in best_chr) print id, best_chr[id], best_pos[id]
        }
    ' \
    | sort --parallel="$THREADS" -t $'\t' -k1,1 > "$TMPDIR/rsid_to_hg19.sorted.tsv"

echo "[$(date)] Preparing input data"
read_input | awk 'BEGIN{FS=OFS="\t"} NR==1{print > "'$TMPDIR'/header.tsv"; next} {print}' \
    > "$TMPDIR/input_data.tsv"

echo "[$(date)] Annotating input with hg19 coordinates (awk left join)"
awk 'BEGIN{FS=OFS="\t"}
    NR==FNR { coord[$1]=$2 FS $3; next }
    {
        rsid=$1
        if (rsid in coord) print coord[rsid], $0
        else print "NA", "NA", $0
    }
' "$TMPDIR/rsid_to_hg19.sorted.tsv" "$TMPDIR/input_data.tsv" > "$TMPDIR/with_coords.unsorted.tsv"

echo "[$(date)] Sorting output by CHR and POS"
{
    awk 'BEGIN{FS=OFS="\t"} {print "CHR","GENPOS",$0}' "$TMPDIR/header.tsv"
    awk 'BEGIN{FS=OFS="\t"}
        {
            c=toupper($1)
            if (c=="X") k=23;
            else if (c=="Y") k=24;
            else if (c=="M" || c=="MT") k=25;
            else if (c ~ /^[0-9]+$/) k=c+0;
            else k=999;
            p=($2 ~ /^[0-9]+$/)?$2+0:999999999;
            print k,p,$0
        }' "$TMPDIR/with_coords.unsorted.tsv" \
        | sort --parallel="$THREADS" -t $'\t' -k1,1n -k2,2n \
        | cut -f3-
} > "$TMPDIR/final.tsv"

if [[ "$OUTFILE" == *.gz ]]; then
    if command -v bgzip >/dev/null 2>&1; then
        bgzip -@ "$THREADS" -c "$TMPDIR/final.tsv" > "$OUTFILE"
    else
        gzip -c "$TMPDIR/final.tsv" > "$OUTFILE"
    fi
else
    cp "$TMPDIR/final.tsv" "$OUTFILE"
fi

echo "[$(date)] Wrote output: $OUTFILE"


# report job stats
sacct -j $SLURM_JOB_ID --format=JobID,JobName,Partition,NodeList,Start,End,Elapsed,State,ExitCode,MaxRSS,MaxVMSize,AveCPU,AveRSS,AveVMSize