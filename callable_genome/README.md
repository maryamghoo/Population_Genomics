# Callable Genome and SNP Filtering Pipeline

## Overview
This workflow identifies callable genomic regions across multiple samples and filters a SNP VCF to retain variants located within commonly callable regions. It consists of four scripts that should be run in order:

1. `1_callable_beds.sh`
2. `2_multiinter_callable.sh`
3. `3_common_callable.sh`
4. `4_callable_vcf.sh`

The workflow starts from BAM files, creates callable regions for each sample, identifies regions callable across a specified number of individuals, and filters the SNP VCF to those regions.

## Purpose
This workflow is used to generate a SNP dataset restricted to reliably callable genomic regions for downstream analyses such as ROH.

The workflow performs:

- Identify callable regions for each sample based on sequencing depth
- Exclude sex chromosomes
- Generate one callable BED file per sample
- Compare callable regions across all samples
- Keep regions callable in at least a specified number of individuals
- Merge common callable regions
- Calculate the total common callable genome size
- Filter a SNP VCF to the common callable genome
- Index the final filtered VCF


## Step 1: Create Callable BED Files

See: `1_callable_beds.sh`

This script reads BAM files from a BAM list and identifies genomic regions with sequencing depth between the specified minimum and maximum depth. Sex chromosomes are excluded. It generates one callable BED file per sample and reports the callable genome size for each sample.

### Expected Inputs
- BAM list containing full paths to BAM files:
```
/data/sample1.bam
/data/sample2.bam
/data/sample3.bam
```
- BAM files
- Sex chromosome list:
```
chrZ
chrW
```

### Output
```
sample1.callable.bed
sample2.callable.bed
sample3.callable.bed
callable_sizes.tsv
```

### Usage
- Update the following in the script:
```
MAXJOBS=4
MIN_DEPTH=2
MAX_DEPTH=21
BAMLIST="/path/to/bamlist.txt"
SEX_LIST="/path/to/sex_chromosomes.txt"
OUTTSV="/path/to/callable_sizes.tsv"
BEDDIR="/path/to/callable_beds"
```
- Submit the job:
`sbatch 1_callable_beds.sh`


## Step 2: Compare Callable Regions Across Samples

See: `2_multiinter_callable.sh`

This script uses `bedtools multiinter` to compare the callable BED files generated in Step 1 and reports how many individuals have each genomic interval as callable.

### Expected Input
A BED list containing paths to the callable BED files:
```
/path/to/sample1.callable.bed
/path/to/sample2.callable.bed
/path/to/sample3.callable.bed
```

### Output
`callable_multiinter.bed`

### Usage
- Update the following in the script:
```
BEDLIST="/path/to/bedlist.txt"
OUT="/path/to/callable_multiinter.bed"
```
- Submit the job:
`sbatch 2_multiinter_callable.sh`


## Step 3: Create Common Callable Regions

See: `3_common_callable.sh`

This script uses the multi-intersection output from Step 2 and retains regions callable in at least the specified number of individuals. These regions are sorted and merged to create the final common callable genome. It also reports basic statistics and the total callable genome size.

### Expected Input
`callable_multiinter.bed`

### Output
`callable_common.bed`
`callable_common.merged.bed`

### Usage
- Update the following in the script:
```
INPUT="/path/to/callable_multiinter.bed"
COMMON="/path/to/callable_common.bed"
MERGED="/path/to/callable_common.merged.bed"
THRESHOLD=218
```
`THRESHOLD` specifies the minimum number of individuals in which a region must be callable.
- Submit the job:
`sbatch 3_common_callable.sh`


## Step 4: Filter SNP VCF to Callable Regions

See: `4_callable_vcf.sh`

This script filters an input SNP VCF using the common callable regions generated in Step 3. Only SNPs located within the common callable genome are retained. The resulting compressed VCF is then indexed.

### Expected Inputs
- SNP VCF:
`input_snps.vcf.gz`
- Common callable BED:
`callable_common.merged.bed`

### Output
`snps_callable.vcf.gz`
`snps_callable.vcf.gz.csi`

### Usage
- Update the following in the script:
```
VCF="/path/to/input_snps.vcf.gz"
BED="/path/to/callable_common.merged.bed"
OUT="/path/to/snps_callable.vcf.gz"
```
- Submit the job:
`sbatch 4_callable_vcf.sh`


## General Usage Notes
- This workflow was run using `bedtools v2.31.1`, `bcftools v1.20`, and `samtools v1.10`. Ensure they are installed and available in your environment.
- Run the four scripts **in order**.
- Adjust the SLURM settings (partition, CPUs, memory, runtime, log paths, etc.) based on the cluster configuration and dataset size.
- Ensure `bedtools` and `bcftools` are available in the environment.
- Input BAM files should be coordinate sorted.
- Update all `/path/to/...` variables before submitting the jobs.
- Check the output and error logs after each step before continuing to the next step.
