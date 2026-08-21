#!/bin/bash -l
#SBATCH -J admixture
#SBATCH -o /admixture_results/logs/admixture.out
#SBATCH -e /admixture_results/logs/admixture.err
#SBATCH -p normal
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=4000
#SBATCH -t 48:00:00

cd /path/to/output_directory

# Convert VCF to PLINK BED format
/path/to/plink \
  --vcf /path/to/snps_afterLD.vcf \
  --allow-extra-chr \
  --make-bed \
  --out adinput


# Run ADMIXTURE for multiple K values
for K in {1..10}; do
  /path/to/admixture \
    -B \
    --cv=10 \
    -j10 \
    adinput.bed $K | tee log${K}.out &

  # Run up to 4 ADMIXTURE jobs simultaneously
  [[ $((K % 4)) -eq 0 ]] && wait
done

wait
