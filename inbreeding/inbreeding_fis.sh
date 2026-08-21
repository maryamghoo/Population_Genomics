# Estimate individual-level observed and expected homozygosity
# and the inbreeding coefficient (F) within each genetic cluster
# using autosomal biallelic SNPs after LD pruning

vcftools \
  --vcf /path/to/snps_afterLD_unrelated.vcf \
  --keep /path/to/<cluster_name>.samples.txt \
  --het \
  --out /path/to/output/<cluster_name>
  