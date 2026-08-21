##### Workflow:
# 1. Calculate IBD within sampling sites
# 2. Identify related pairs using PI_HAT
# 3. Remove one individual from each related pair
# 4. Generate unrelated SNP datasets before and after LD pruning

##### USAGE:
# Update all file paths, filenames, output names, and SNP ID lists according to your project structure.


##### IMPORTANT:
# The .irem file used with --keep should contain sample IDs from a specific sampling site.

##### Example format:

# sample_01 sample_01
# sample_02 sample_02
# sample_03 sample_03

# PLINK expects two columns:
# Family ID (FID) and Individual ID (IID).
# If family IDs are unavailable, the sample ID is commonly repeated in both columns.


##### evaluate relatedness within each sampling site #####

# calculate identity-by-descent (IBD) statistics for individuals within a specific sampling site

/path/to/plink \
  --vcf /path/to/snps_afterLD.vcf \
  --allow-extra-chr \
  --keep /path/to/<site_name>.samples.irem \
  --genome \
  --out /path/to/output/relatedness_<site_name>


##### main output file
# relatedness_<site_name>.genome

# the .genome file contains pairwise relatedness statistics,
# including the PI_HAT value for each pair of individuals

# individuals with PI_HAT >= 0.5 were considered closely related.
# for each related pair, one individual was randomly removed.

# manually create a text file listing the IDs of individuals selected for removal:
# ibd_rmv_samples.txt

# example format:
# sample_01
# sample_15
# sample_22


##### remove related individuals from the VCF file #####

# remove individuals listed in ibd_rmv_samples.txt and create a new vcf file

bcftools view \
  -S ^/path/to/ibd_rmv_samples.txt \
  /path/to/passedvariants.vcf \
  -o /path/to/unrelated_individuals.vcf \
  -O v


##### generate SNP ID lists #####

# The SNP ID lists generated below are used to reconstruct
# VCF files containing only unrelated individuals while
# preserving the SNP sets before and after LD pruning.

# SNPs before LD pruning
# extract SNP IDs from the .bim file (2nd column)

cut -f2 /path/to/snps_beforeLD.bim > \
/path/to/snps_beforeLD.txt

# SNPs retained after LD pruning
# PLINK generates snps_afterLD.prune.in automatically

cp /path/to/snps_afterLD.prune.in \
/path/to/snps_afterLD.txt


##### generate SNP dataset before LD pruning #####

# keep SNPs before LD pruning
# using SNP IDs from a previous filtering step (extracted from the PLINK .bim file)

bcftools view \
  -i 'ID=@/path/to/snps_beforeLD.txt' \
  /path/to/unrelated_individuals.vcf \
  -o /path/to/snps_beforeLD_unrelated.vcf \
  -O v
  
  
##### generate SNP dataset after LD pruning #####

# keep only SNPs retained after LD pruning
# using SNP IDs from a previous filtering step

bcftools view \
  -i 'ID=@/path/to/snps_afterLD.txt' \
  /path/to/unrelated_individuals.vcf \
  -o /path/to/snps_afterLD_unrelated.vcf \
  -O v
