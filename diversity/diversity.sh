# Estimate nucleotide diversity (pi) for each sampling region
# using autosomal biallelic SNPs before LD pruning

vcftools \
  --vcf /path/to/snps_beforeLD_unrelated.vcf \
  --keep /path/to/<site_name>.samples.txt \
  --window-pi 50000 \
  --window-pi-step 50000 \
  --out /path/to/output/pi_<site_name>


# Estimate observed and expected heterozygosity for each sampling region
# using autosomal biallelic SNPs after LD pruning

vcftools \
  --vcf /path/to/snps_afterLD_unrelated.vcf \
  --keep /path/to/<site_name>.samples.txt \
  --het \
  --out /path/to/output/het_<site_name>
  