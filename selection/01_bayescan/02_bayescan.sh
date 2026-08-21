#!/bin/bash
#SBATCH -J bayescan
#SBATCH -o /bayescan/logs/bayescan.out
#SBATCH -e /bayescan/logs/bayescan.err
#SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH -t 20-00:00:00
#SBATCH --mem=64G

# Paths
INDIR="/path/to/bayescan"
OUTDIR="/path/to/bayescan/bayescan_output"

mkdir -p "$OUTDIR"

# Run BayeScan
/path/to/bayescan_2.1 "$INDIR/bayescan_input.txt" \
    -od "$OUTDIR" \
    -snp \
    -pr_odds 10000 \
    -threads "$SLURM_CPUS_PER_TASK" \
    -n 5000 \
    -thin 20 \
    -nbp 10 \
    -pilot 2000 \
    -burn 20000
