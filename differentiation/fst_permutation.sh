#!/bin/bash -l
#SBATCH -J fst_permutation
#SBATCH -o /path/to/logs/fst_permutation.out
#SBATCH -e /path/to/logs/fst_permutation.err
#SBATCH -p normal
#SBATCH --cpus-per-task=1
#SBATCH -t 48:00:00

# Stop on errors, undefined variables, and failed pipes
set -euo pipefail

# Define population/sampling-site names
POP1_NAME="<site_1>"
POP2_NAME="<site_2>"

# Define input files and directories
VCF="/path/to/snps_beforeLD_unrelated.vcf"
SAMPLE_DIR="/path/to/sample_lists"

# Define output directories
OUTDIR="/path/to/fst_results"
PERM_DIR="${OUTDIR}/permutations"

# Create output directories if they do not exist
mkdir -p "${OUTDIR}"
mkdir -p "${PERM_DIR}"

# Number of permutations
N_PERM=100

# Number of parallel jobs
N_JOBS=8

# 1. Calculate observed FST 

# Calculate windowed FST between the two original sampling sites
vcftools \
  --vcf "${VCF}" \
  --weir-fst-pop "${SAMPLE_DIR}/${POP1_NAME}.txt" \
  --weir-fst-pop "${SAMPLE_DIR}/${POP2_NAME}.txt" \
  --fst-window-size 50000 \
  --fst-window-step 25000 \
  --out "${OUTDIR}/${POP1_NAME}.${POP2_NAME}"


# 2. Get original population sizes

POP1_COUNT=$(wc -l < "${SAMPLE_DIR}/${POP1_NAME}.txt")
POP2_COUNT=$(wc -l < "${SAMPLE_DIR}/${POP2_NAME}.txt")


# 3. Create CSV file for permutation results

echo "permutation,average_weighted_fst,average_mean_fst" \
  > "${PERM_DIR}/${POP1_NAME}.${POP2_NAME}.average_fst_permutation.csv"


# 4. Run FST permutations

for i in {1..100}
do

  # Randomly shuffle individuals between the two sampling sites
  cat \
    "${SAMPLE_DIR}/${POP1_NAME}.txt" \
    "${SAMPLE_DIR}/${POP2_NAME}.txt" | \
    shuf > "${PERM_DIR}/${POP1_NAME}.${POP2_NAME}.shuffled.txt"

  # Split shuffled individuals while preserving original sample sizes
  head -n "${POP1_COUNT}" \
    "${PERM_DIR}/${POP1_NAME}.${POP2_NAME}.shuffled.txt" \
    > "${PERM_DIR}/perm_${POP1_NAME}_${POP2_NAME}.txt"

  tail -n "${POP2_COUNT}" \
    "${PERM_DIR}/${POP1_NAME}.${POP2_NAME}.shuffled.txt" \
    > "${PERM_DIR}/perm_${POP2_NAME}_${POP1_NAME}.txt"


  # Calculate FST for permuted sampling-site assignments

  VCFTOOLS_OUTPUT=$(vcftools \
    --vcf "${VCF}" \
    --weir-fst-pop "${PERM_DIR}/perm_${POP1_NAME}_${POP2_NAME}.txt" \
    --weir-fst-pop "${PERM_DIR}/perm_${POP2_NAME}_${POP1_NAME}.txt" \
    --fst-window-size 50000 \
    --fst-window-step 25000 \
    --stdout 2>/dev/null)


  # Calculate average permutation FST and append to CSV
  
  echo "${VCFTOOLS_OUTPUT}" | \
  	awk -F"\t" \
  	'BEGIN {sum5=0; sum6=0; n=0}
  	NR>1 && $5 != "nan" && $6 != "nan" {
  	sum5 += $5;
     sum6 += $6;
     n++
   	}
   END {
     if (n > 0)
       print '"$i"' "," sum5/n "," sum6/n;
     else
       print '"$i"' ",NA,NA"
   }' \
  >> "${PERM_DIR}/${POP1_NAME}.${POP2_NAME}.average_fst_permutation.csv"

done
