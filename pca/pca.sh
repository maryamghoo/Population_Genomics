#!/bin/bash

# all regions (no subset)
/path/to/plink \
  --vcf /path/to/snps_afterLD.vcf \
  --allow-extra-chr \
  --pca 20 \
  --out /path/to/output/all_regions

# Convert eigenvalues to CSV: PC1,eigenvalue ... PC20,eigenvalue
awk '{printf "PC%02d,%s\n", NR, $1}' \
  /path/to/output/all_regions.eigenval \
  > /path/to/output/all_regions.eigenvalues.csv

awk '{$1=$1}1' OFS=',' \
  /path/to/output/all_regions.eigenvec \
  > /path/to/output/all_regions.eigenvectors.csv

# site-specific PCA (mainland_bc, cbc_sbc, etc.)
/path/to/plink \
  --vcf /path/to/snps_afterLD.vcf \
  --keep /path/to/<site>_samples.txt \
  --allow-extra-chr \
  --pca 20 \
  --out /path/to/output/<site>

awk '{printf "PC%02d,%s\n", NR, $1}' \
  /path/to/output/<site>.eigenval \
  > /path/to/output/<site>.eigenvalues.csv

awk '{$1=$1}1' OFS=',' \
  /path/to/output/<site>.eigenvec \
  > /path/to/output/<site>.eigenvectors.csv
  