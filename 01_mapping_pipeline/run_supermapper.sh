#!/bin/bash -l
#SBATCH -J SAMPLE_ID_mapping
#SBATCH -o /mapping_results/logs/SAMPLE_ID.out
#SBATCH -e /mapping_results/logs/SAMPLE_ID.err
#SBATCH -p normal
#SBATCH --cpus-per-task=12
#SBATCH -t 48:00:00
#SBATCH --mem-per-cpu=8000

cd /path/to/output_directory

# Define sample name
SAMPLE_NAME="SAMPLE_ID"

# Get unmapped read using BWA aligner
supermapper \
-i /data/fastq/${SAMPLE_NAME}_R1.fastq.gz \
-j /data/fastq/${SAMPLE_NAME}_R2.fastq.gz \
-r /reference/reference_genome.fna \
-t 12 \
-o ${SAMPLE_NAME} \
-g "@RG\tID:$SAMPLE_NAME\tSM:$SAMPLE_NAME\tPL:ILLUMINA"
