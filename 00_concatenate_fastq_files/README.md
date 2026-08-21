# FASTQ Files Concatenation Script

## Overview
This repository contains the script concatenateing FASTQ files generated from multiple sequencing lanes into a single pair of files per sample (`R1` and `R2`). It is used as a preprocessing step before running downstream pipelines.
See: `concatenate_fastq_files.sh`

## Purpose
Many sequencing runs split reads across multiple lanes (e.g., L001–L004). This script merges those lane-specific FASTQ files into unified files so they can be processed together.

The output of this step is used as input for the `run_supermapper.sh` pipeline, which performs:
- Read trimming
- Quality cleaning
- Alignment

## Expected inputs
The script expects FASTQ files in the following format:
```
SAMPLE_ID_L001_R1_001.fastq.gz
SAMPLE_ID_L002_R1_001.fastq.gz
SAMPLE_ID_L003_R1_001.fastq.gz
SAMPLE_ID_L004_R1_001.fastq.gz

SAMPLE_ID_L001_R2_001.fastq.gz
SAMPLE_ID_L002_R2_001.fastq.gz
SAMPLE_ID_L003_R2_001.fastq.gz
SAMPLE_ID_L004_R2_001.fastq.gz
```

## Output
For each sample, the script generates:
```
SAMPLE_ID_R1.fastq.gz
SAMPLE_ID_R2.fastq.gz
```

## Usage
- Adjust the SLURM job settings (partition, CPUs, memory, runtime, output paths, etc.) based on your cluster configuration, dataset size, software requirements, and expected workload.
- Update the following in the script:
```
concat_results/logs/SAMPLE_ID.out
concat_results/logs/SAMPLE_ID.err
input_dir=/path/to/fastq_files
output_dir=/path/to/output_directory
SAMPLE_NAME="SAMPLE_ID"
```  
- Run the script once per sample (update all SAMPLE_IDs each time):
`sbatch concatenate_fastq_files.sh`
