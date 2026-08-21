# Downsampling BAM Files Script

## Overview
This repository contains the script for downsampling BAM alignment files to a maximum number of reads using samtools, then generates indexed, sorted BAM files and alignment statistics.
See: `downsampling.sh`

## Purpose
Sequencing samples often have uneven read depth across individuals. This script standardizes BAM files by randomly subsampling reads to a defined maximum read count, helping reduce coverage bias in downstream analyses.

The script performs:

- Estimate total mapped reads
- Calculate downsampling fraction
- Randomly subsample reads
- Generate BAM statistics
- Sort BAM files
- Index output BAM files

## Expected Inputs
- The script expects a text file containing full paths to BAM files:
```
/path/to/sample1.bam
/path/to/sample2.bam
/path/to/sample3.bam
```
- Each BAM file should already have an index file:
`sample.bam.bai`

## Output
For each sample, the script generates:
```
sample_subsampled.bam
sample_subsampled.bam.bai
sample_subsampled.bamstats
sample_subsampled.sorted.bam
sample_subsampled.sorted.bam.bai
```

## Usage
- This script was run using `samtools v1.10`. Ensure `samtools` is installed and available in your environment.
- Adjust the SLURM job settings (partition, CPUs, memory, runtime, output paths, etc.) based on your cluster configuration, dataset size, software requirements, and expected workload.
- Update the following in the script:
```
/downsampling_results/downsampling.out
/downsampling_results/downsampling.err
max_reads=250000000
input_list=/path/to/downsampling_target_bamlist.txt
output_dir=/path/to/downsampled/directory
```
- Submit the job:
`sbatch downsampling.sh`
