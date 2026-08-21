#!/bin/bash -l
#SBATCH -J variantCalling
#SBATCH -o /variant_calling_results/logs/variants.out
#SBATCH -e /variant_calling_results/logs/variants.err
#SBATCH -p normal
#SBATCH --cpus-per-task=64
#SBATCH -t 48:00:00

cd /path/to/output_directory

# Calling Variants
/path/to/bcftools \
mpileup \
-C 50 \
-q 20 \
-Q 20 \
--threads 64 \
-Ou \
--per-sample-mF \
--annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \
-f /reference/reference_genome.fna \
-b /data/bamlist.txt | \
/path/to/bcftools \
call \
--threads 64 \
-A \
-f GQ \
-mv -Oz > variants.vcf.gz

# Index VCF
/path/to/bcftools index \
--threads 16 \
variants.vcf.gz

# Unzip VCF without overwriting original gz file
gunzip -c variants.vcf.gz > variants.vcf

# Generate stats file
/path/to/bcftools stats \
variants.vcf.gz > variants.stats
