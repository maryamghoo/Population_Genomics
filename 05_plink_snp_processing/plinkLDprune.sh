#!/bin/bash -l
#SBATCH -J plinkLDprune
#SBATCH -o /plink_results/logs/plink_ld.out
#SBATCH -e /plink_results/logs/plink_ld.err
#SBATCH -p normal
#SBATCH --cpus-per-task=16
#SBATCH -t 48:00:00

cd /path/to/output_directory

# SNP filtering + LD pruning list generation
/path/to/plink \
--vcf /path/to/passedvariants.vcf \
--allow-extra-chr \
--vcf-filter \
--snps-only \
--biallelic-only strict list \
--mind 0.2 \
--geno 0.2 \
--maf 0.02 \
--make-bed \
--not-chr JBGGOT010000034.1 JBGGOT010000035.1 JBGGOT010000036.1 \
JBGGOT010000037.1 JBGGOT010000038.1 JBGGOT010000039.1 JBGGOT010000040.1 \
JBGGOT010000041.1 JBGGOT010000042.1 \
--indep-pairwise 50 kb 5 0.5 \
--out snps_afterLD

# Keep only pruned SNPs and export to VCF
/path/to/plink \
--bfile snps_afterLD \
--allow-extra-chr \
--extract snps_afterLD.prune.in \
--recode vcf \
--out snps_afterLD

# Generate stats
bcftools stats snps_afterLD.vcf > snps_afterLD.stats
