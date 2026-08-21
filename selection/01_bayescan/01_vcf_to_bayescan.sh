#!/bin/bash
#SBATCH -J vcf2bayescan
#SBATCH -o /vcf2bayescan/logs/vcf2bayescan.out
#SBATCH -e /vcf2bayescan/logs/vcf2bayescan.err
#SBATCH -p normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH -t 5-00:00:00

# Paths
VCF="/path/to/snps_afterLD_unrelated.vcf"
POP_FILE="/path/to/sample_populations.txt"

OUTFILE="/path/to/bayescan_input.txt"

perl /path/to/vcf2bayescan.pl \
    -infile "$VCF" \
    -popFile "$POP_FILE" \
    -colNum 1 \
    -outfile "$OUTFILE"
