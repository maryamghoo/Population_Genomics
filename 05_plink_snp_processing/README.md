# PLINK SNP Filtering and LD Pruning Scripts

## Overview
These scripts process a filtered VCF file using PLINK to generate high-quality biallelic autosomal SNP dataset, remove low-quality samples and loci, apply missingness and minor allele frequency filters, and optionally perform linkage disequilibrium (LD) pruning. Final datasets are exported as VCF files with summary statistics generated using bcftools.
See:
`plink_snp_filter.sh`
`plink_ld_prune.sh`

## Purpose
After raw variant calling and variant filtering, additional SNP-level filtering is commonly required before downstream analyses such as:

Relatedness analysis
Population structure
Diversity estimates

These scripts create:

- SNP dataset before LD pruning
- LD-pruned SNP dataset

## Expected Inputs
Both scripts expect:
A filtered PASS variant file from the previous pipeline step. 
`passedvariants.vcf`

## Script 1: SNP Filtering Before LD Pruning

This script performs:
- Keep SNPs only
- Keep strictly biallelic loci only
- Remove individuals with >20% missing genotypes (--mind 0.2)
- Remove SNPs with >20% missingness (--geno 0.2)
- Remove rare SNPs (--maf 0.02)
- Remove excluded sex chromosomes/scaffolds
- Export filtered VCF
- Generate summary statistics

### Output
`snps_beforeLD.bed`
`snps_beforeLD.bim`
`snps_beforeLD.fam`
`snps_beforeLD.vcf`
`snps_beforeLD.stats`
`snps_beforeLD.log`
`snps_beforeLD.irem`
`snps_beforeLD.nosex` generated when sex information is missing or undefined
`snps_beforeLD-temporary.skip.3allele` listing skipped multiallelic variants

## Script 2: SNP Filtering + LD Pruning

This script performs all filters above, then additionally removes correlated SNPs using a 50 kb window, step size of 5 SNPs, and LD threshold of r² = 0.5.

### Output
`snps_afterLD.bed`
`snps_afterLD.bim`
`snps_afterLD.fam`
`snps_afterLD.prune.in`
`snps_afterLD.prune.out`
`snps_afterLD.vcf`
`snps_afterLD.stats`
`snps_afterLD.log`
`snps_afterLD.irem`
`snps_afterLD.nosex` generated when sex information is missing or undefined
`snps_afterLD-temporary.skip.3allele` listing skipped multiallelic variants

## Usage
- These script were run using `plink v1.9` and `bcftools v1.20`. Ensure `plink` and `bcftools` are installed and available in your environment.
- Adjust the SLURM job settings (partition, CPUs, memory, runtime, output paths, etc.) based on your cluster configuration, dataset size, software requirements, and expected workload.
- Adjust SNP filtering thresholds, LD pruning parameters, and the list of excluded chromosomes/scaffolds based on your species, reference genome assembly, marker density, sample size, sequencing quality, and downstream analysis goals.
- Update the following in each script:
```
/plink_results/logs/plink_snps.out
/plink_results/logs/plink_snps.err
/plink_results/logs/plink_ld.out
/plink_results/logs/plink_ld.err
/path/to/output_directory
/path/to/plink
/path/to/passedvariants.vcf
```
- Submit jobs:
`sbatch plink_snp_filter.sh`
`sbatch plink_ld_prune.sh`
