# Variant Filtering Script

## Overview
This script filters raw variant calls from a VCF file using bcftools. It applies multiple quality, depth, strand-bias, and support-read filters to retain higher-confidence variants, extracts only PASS variants, renames variant IDs, and generates summary statistics.
See: `variant_filtering.sh`

## Purpose
Raw variant calls often contain low-quality or weakly supported sites. This script removes or flags unreliable variants to generate a cleaner dataset for downstream population genetics, association, or diversity analyses.

The script performs:

- Read input VCF file
- Apply quality score filtering
- Apply mapping quality filtering
- Filter extreme depth values
- Remove strand/read-position bias signals
- Require minimum alternate allele support reads
- Mark failing variantswith `FAIL`
- Keep only `PASS` variants
- Rename variant IDs sequentially (`SNP_1`, `SNP_2`,...)
- Generate summary statistics for final variants
- Output filtered VCF file

## Expected Inputs
The script expects:
- Input VCF file generated from a previous variant calling step.
`variants.vcf`

## Output
The script generates:
`filteredvariants.vcf`
`passedvariants.vcf`
`passedvariants.stats`

## Usage
- This script was run using `bcftools v1.20`. Ensure `bcftools` is installed and available in your environment.
- Adjust the SLURM job settings (partition, CPUs, memory, runtime, output paths, etc.) based on your cluster configuration, dataset size, software requirements, and expected workload.
- The filtering thresholds may need adjustment depending on species, sequencing depth, and study design.
- Update the following in the script:
```
/variant_filtering_results/logs/filtered_variants.out
/variant_filtering_results/logs/filtered_variants.err
/path/to/output_directory
/path/to/bcftools
/path/to/variants.vcf
```
- Submit the job:
`sbatch variant_filtering.sh`
