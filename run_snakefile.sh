#!/bin/bash
#SBATCH --job-name=RNA_seq
#SBATCH --output=logs/RNA_seq.out
#SBATCH --error=logs/RNA_seq.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=120:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

# Exit on error
set -e

# Create necessary directories
mkdir -p logs/cluster_logs
mkdir -p results/fastqc
mkdir -p results/trimmed
mkdir -p results/star
mkdir -p results/counts
mkdir -p results/multiqc

# Set working directory
cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Federica_RNA

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Unlock snakemake working directory if necessary
# snakemake --unlock

# Run snakemake with dry-run first to validate
echo "Performing dry-run to validate workflow..."
snakemake --snakefile Snakefile --dry-run --rerun-incomplete

if [ $? -eq 0 ]; then
    echo "Dry-run successful, starting actual run..."
    # Run snakemake
    snakemake \
        --snakefile Snakefile \
        --executor slurm \
        --rerun-incomplete \
        --jobs 100 \
        --default-resources \
            slurm_partition=workq \
            mem_mb=32000 \
            runtime=1440 \
            threads=8 \
            nodes=1 \
        --jobscript slurm-jobscript.sh \
        --latency-wait 60 \
        --rerun-incomplete \
        --keep-going
else
    echo "Dry-run failed, please check your Snakefile for errors"
    exit 1
fi
