#!/bin/bash -l
#SBATCH -J nonsyn_freq
#SBATCH -o /nonsyn_freq/logs/nonsyn_freq.out
#SBATCH -e /nonsyn_freq/logs/nonsyn_freq.err
#SBATCH -p normal
#SBATCH --cpus-per-task=1
#SBATCH -t 24:00:00

set -euo pipefail

##### Nonsynonymous SNP allele-frequency analysis

# Define directories
DATA_DIR="bash_data"
OUTDIR="bash_results"

mkdir -p "${OUTDIR}"
mkdir -p "${OUTDIR}/frequencies"

# Define input files
ANNOTATION="${DATA_DIR}/fdr_outliers_annotated.fsj_multianno.txt"
VCF="/path/to/snps_afterLD_unrelated.vcf.gz"
SAMPLE_DIR="${DATA_DIR}/sample_lists"

# Define genetic clusters

POPULATIONS=(
  "archbold"
  "cbc"
  "sbc"
  "sl"
  "ocala"
  "merritt"
)

# 1. Identify nonsynonymous SNPs from ANNOVAR results
awk -F'\t' \
  'BEGIN{OFS="\t"} NR > 1 && $9 == "nonsynonymous SNV" {print $1,$2,$3}' \
  "${ANNOTATION}" \
  > "${OUTDIR}/nonsynonymous_snps.txt"

# 2. Extract nonsynonymous SNPs from the VCF
bcftools view \
  -R "${OUTDIR}/nonsynonymous_snps.txt" \
  "${VCF}" \
  -Oz \
  -o "${OUTDIR}/nonsynonymous_snps.vcf.gz"

bcftools index \
  -f \
  "${OUTDIR}/nonsynonymous_snps.vcf.gz"

# 3. Calculate allele frequencies for each genetic cluster
for POP in "${POPULATIONS[@]}"
do

  vcftools \
    --gzvcf "${OUTDIR}/nonsynonymous_snps.vcf.gz" \
    --keep "${SAMPLE_DIR}/${POP}.txt" \
    --freq \
    --out "${OUTDIR}/frequencies/${POP}_freq"

done

# 4. Combine allele-frequency results
COMBINED="${OUTDIR}/all_pop_freq.txt"

echo -e "cluster\tCHROM\tPOS\tN_ALLELES\tN_CHR\tALLELE1\tALLELE2" \
  > "${COMBINED}"

for POP in "${POPULATIONS[@]}"
do

  awk -v pop="${POP}" \
    'BEGIN{OFS="\t"} NR > 1 {print pop,$1,$2,$3,$4,$5,$6}' \
    "${OUTDIR}/frequencies/${POP}_freq.frq" \
    >> "${COMBINED}"

done
