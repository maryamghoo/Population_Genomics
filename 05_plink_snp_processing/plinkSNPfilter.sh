#!/bin/bash -l
#SBATCH -J plinkSNPfilter
#SBATCH -o /plink_results/logs/plink_snps.out
#SBATCH -e /plink_results/logs/plink_snps.err
#SBATCH -p normal
#SBATCH --cpus-per-task=16
#SBATCH -t 48:00:00

cd /path/to/output_directory

# Create BED/BIM/FAM dataset with SNP filtering
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
--out snps_beforeLD

# Export filtered data to VCF
/path/to/plink \
--bfile snps_beforeLD \
--allow-extra-chr \
--recode vcf \
--out snps_beforeLD

# Generate stats
bcftools stats snps_beforeLD.vcf > snps_beforeLD.stats
