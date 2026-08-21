# Mapping with Supermapper 

## Overview
This repository contains the mapping script used for paired-end whole genome sequencing reads from the Florida Scrub-Jay (Aphelocoma coerulescens) in our study.
Reads were processed using `supermapper`, a lab-developed wrapper pipeline that performs read cleaning/trimming and alignment to a reference genome.
See: `run_supermapper.sh`

## Important note
`supermapper` is a local lab script developed by Dr. Robert Fitak, Department of Biology; Genomics and Bioinformatics Cluster, University of Central Florida
robert.fitak@ucf.edu
https://github.com/rfitak
   
Users should either:
1. obtain access to the original `supermapper` script from the lab, or
2. reproduce the same workflow with equivalent tools for trimming, alignment, and BAM generation.

In our case, `supermapper` was run on a Linux cluster through a SLURM job script.

## Expected inputs
- Paired-end FASTQ files (`R1` and `R2`)
- Reference genome FASTA file
The script expects FASTQ files in the following format:
```
SAMPLE_ID_R1.fastq.gz
SAMPLE_ID_R2.fastq.gz
```

## Output
For each sample, the script generates:
```
SAMPLE_ID.bamstats
SAMPLE_ID.clean.sorted.bam
SAMPLE_ID.clean.sorted.bam.bai
SAMPLE_ID.html
```

## Usage
- Adjust the SLURM job settings (partition, CPUs, memory, runtime, output paths, etc.) based on your cluster configuration, dataset size, software requirements, and expected workload.
- Update the following in the script:
```
mapping_results/logs/SAMPLE_ID.out
mapping_results/logs/SAMPLE_ID.err
/path/to/output_directory
data/fastq/${SAMPLE_NAME}_R1.fastq.gz
data/fastq/${SAMPLE_NAME}_R2.fastq.gz
reference/reference_genome.fna
${SAMPLE_NAME}
SAMPLE_NAME="SAMPLE_ID"
```  
- Run the script separately for each sample (update SAMPLE_IDs each time):
`sbatch run_supermapper.sh`
