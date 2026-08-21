# Relatedness Filtering Workflow

## Overview
This script documents the workflow used to identify and remove closely related individuals from a SNP dataset. It combines commands from PLINK and bcftools to calculate pairwise relatedness, remove selected individuals from a VCF file, and generate VCF files containing only unrelated individuals before and after LD pruning.
See: `remove_related_individuals.sh`

## Purpose

Closely related individuals can bias many downstream population genomic analyses. This workflow identifies related pairs based on identity-by-descent (IBD) estimates, removes one individual from each highly related pair, and reconstructs SNP datasets containing only unrelated individuals.

The workflow performs:

- Calculate pairwise IBD statistics within each sampling site
- Identify closely related individuals using PI_HAT values
- Remove one individual from each related pair
- Remove selected individuals from the filtered VCF file
- Generate SNP ID lists before and after LD pruning
- Create VCF files containing only unrelated individuals for both SNP datasets

## Expected Inputs
The workflow expects:
- Filtered VCF generated from the previous variant filtering step.
`passedvariants.vcf`
- Sampling-site sample list(s).
`<site_name>.samples.irem`
- PLINK SNP datasets generated from previous SNP filtering steps.
`snps_beforeLD.bim`
`snps_afterLD.prune.in`

## Output
The workflow generates:
`relatedness_<site_name>.genome`
`ibd_rmv_samples.txt`
`unrelated_individuals.vcf`
`snps_beforeLD.txt`
`snps_afterLD.txt`
`snps_beforeLD_unrelated.vcf`
`snps_afterLD_unrelated.vcf`

## Usage
- This workflow was run using `plink v1.9` and `bcftools v1.20`. Ensure `plink` and `bcftools` are installed and available in your environment.
- This file documents a sequence of commands used in the study and is intended as a workflow reference rather than a standalone executable pipeline or SLURM job.
- Follow the commands in order, as later steps depend on files generated or prepared in earlier steps.
- Update all file paths, filenames, sample lists, SNP ID lists, and output names according to your project structure before running the commands.
- After reviewing the `relatedness_<site_name>.genome` file, manually identify pairs with PI_HAT ≥ 0.5, randomly select one individual from each pair for removal, and create `ibd_rmv_samples.txt` before continuing with the remaining steps.
