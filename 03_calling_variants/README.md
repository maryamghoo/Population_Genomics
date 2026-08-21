# Variant Calling Script

## Overview
This script performs SNP/variant calling from multiple BAM alignment files using bcftools. It uses `mpileup` to generate genotype likelihood data and pipes the output directly into `bcftools call` to produce a compressed VCF file, followed by indexing, uncompressing the VCF, and generating summary statistics for downstream analyses.
See: `variant_calling.sh`

## Purpose
This workflow is used after read alignment to detect genomic variants across multiple samples.

The script performs:

- Read aligned BAM files from a BAM list
- Generate pileup/genotype likelihood data
- Apply mapping and base quality filters
- Call SNPs/variants using multiallelic calling mode
- Output compressed VCF results (`.vcf.gz`)
- Index the final VCF file (`.csi` or `.tbi`)
- Create an uncompressed VCF copy (`.vcf`)
- Generate variant statistics report

## Expected Inputs
The script expects:
- BAM list file
- A text file containing full paths to BAM files:
```
/data/sample1.bam
/data/sample2.bam
/data/sample3.bam
```
- Reference genome (should also be indexed)
- Input BAM files should be indexed and present in the BAM files directory.
`sample.bam.bai`

## Output
The script generates a compressed VCF file containing called variants.
`variants.vcf.gz`
`variants.vcf.gz.csi`
`variants.vcf`
`variants.stats`

## Usage
- This script was run using `bcftools v1.20`. Ensure `bcftools` is installed and available in your environment.
- Adjust the SLURM job settings (partition, CPUs, memory, runtime, output paths, etc.) based on your cluster configuration, dataset size, software requirements, and expected workload.
- Update the following in the script:
```
/variant_calling_results/logs/variants.out
/variant_calling_results/logs/variants.err
/path/to/output_directory
/path/to/bcftools
/reference/reference_genome.fna
/data/bamlist.txt
```
- Submit the job:
`sbatch variant_calling.sh`
