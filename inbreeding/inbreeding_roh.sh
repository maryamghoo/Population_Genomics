#!/bin/bash
#SBATCH -J roh
#SBATCH -o /roh/logs/roh.out
#SBATCH -e /roh/logs/roh.err
#SBATCH -p normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -t 3-00:00:00

# Stop on errors, undefined variables, and failed pipes
set -euo pipefail

##### Settings #####

PLINK="/path/to/plink"
VCF="/path/to/snps_callable.vcf.gz"

OUT="/path/to/output/roh"

##########

$PLINK --vcf "$VCF" \
--allow-extra-chr \
--homozyg-window-snp 50 \
--homozyg-kb 1000 \
--homozyg-gap 1000 \
--homozyg-window-het 1 \
--homozyg-window-missing 5 \
--threads 4 \
--out "$OUT"
